import os
import requests
import urllib.request
import itertools
import numpy as np
from pydna.genbank import Genbank
from Bio import SeqIO
from Bio.Restriction import AllEnzymes

from Bio import Entrez
from Bio import Entrez

from Bio.SeqIO import read as seq_read
from io import StringIO


def feature_table(sequence):
    sorted_feat = sorted(sequence.features, key=lambda x: int(x.location.start))
    
    # misc_feature is MCS?
    print(f"{'feature':15} {'start':6}{'end':6} {'notes'}")
    for feature in sorted_feat:
        start = feature.location.start
        end = feature.location.end
        print(f"{feature.type:15} {start:6d}{end:6d} {feature.qualifiers.get("note")}")
        #print(feature.type, feature.location, feature.qualifiers.get("note"))


def load_from_genebank(email, sequence):
    gb = Genbank(users_email=email)  
    seq = gb.nucleotide(sequence)   
    return seq


def cut_enzyme_info(seq, enzyme_list=AllEnzymes):
    cut_sites = {}
    for enzyme in enzyme_list:
        cut = enzyme.search(seq)
        if cut:
            cut_sites[enzyme] = cut
    return cut_sites


def _get_mcs_cuts(sequence, mcs_feat):
    # Get location of MCS
    locs = mcs_feat.location
    start, end = locs.start, locs.end 
    mcs = sequence.seq[start: end]
    
    # Find enzyme cuts
    enz_cuts = cut_enzyme_info(mcs)

    # Get enzymes in the MCS
    mcs_enz = list(enz_cuts.keys())
    
    # Get cuts and add offset
    mcs_cuts = np.array(list(set(itertools.chain.from_iterable(list(enz_cuts.values()))))) + start
    return {"enzymes": mcs_enz, "cuts": mcs_cuts}


def get_mcs_cuts(sequence):
    mcs_list = list(filter(lambda feat: feat.type=='misc_feature', sequence.features))
    mcs_enz_cut = []
    for mcs in mcs_list:
        mcs_enz_cut.append(_get_mcs_cuts(sequence, mcs))
    return mcs_enz_cut



def get_non_mcs_regions(sequence):
    # get non-mcs regions
    mcs_list = list(filter(lambda feat: feat.type=='misc_feature', sequence.features))
    mcs_list
    sequence = sequence
    landmarks = []
    for mcs in mcs_list: 
        locs = mcs.location
        start, end = locs.start, locs.end
        landmarks.append(start)
        landmarks.append(end) 
    
    landmarks = np.array(sorted(landmarks))
    landmarks = np.roll(landmarks, shift=1)
    landmarks
    
    non_mcs_regions = []
    for i in range(0, len(landmarks), 2):
        start, end = landmarks[i:i+2]
        non_mcs = sequence.seq[start: end]
        non_mcs_regions.append(non_mcs)

    return non_mcs_regions



def fetch_species_genome_from_gene_biopython(gene_seq):
    """
    Given a gene Dseqrecord, fetch the complete chromosomal genome of its species using Entrez.
    """
    species_name = gene_seq.annotations["organism"]
    
    # More specific search to get chromosomal genome
    search_terms = [
        f'"{species_name}"[Organism] AND "complete genome"[Title] AND chromosome[Title]',
        f'"{species_name}"[Organism] AND "complete genome"[Title] AND "chromosome 1"[Title]',
        f'"{species_name}"[Organism] AND "complete genome"[Title]',
        f'"{species_name}"[Organism] AND complete genome'
    ]
    
    genome_id = None
    for term in search_terms:
        try:
            handle = Entrez.esearch(db="nucleotide", term=term, retmax=5)
            record = Entrez.read(handle)
            handle.close()
            
            if record["IdList"]:
                # Filter out plasmids and prefer chromosomal entries
                for id in record["IdList"]:
                    summary_handle = Entrez.esummary(db="nucleotide", id=id)
                    summary = Entrez.read(summary_handle)
                    summary_handle.close()
                    
                    title = summary[0]["Title"].lower()
                    # Skip plasmid entries
                    if "plasmid" not in title and "plasmids" not in title:
                        genome_id = id
                        break
                
                if genome_id:
                    break
                    
        except Exception as e:
            print(f"Search failed for term '{term}': {e}")
            continue
    
    if not genome_id:
        raise ValueError(f"Could not find complete chromosomal genome for {species_name}")
    
    # Fetch the genome sequence
    handle = Entrez.efetch(db="nucleotide", id=genome_id, rettype="gb", retmode="text")
    genome_record = seq_read(StringIO(handle.read()), "genbank")
    handle.close()
    
    return genome_record



# Alternative: Use assembly database for more reliable genome fetching
def fetch_genome_from_assembly(gene_seq):
    """
    Fetch complete genome using Assembly database (more reliable)
    """
    species_name = gene_seq.annotations["organism"]
    
    # First search for assemblies
    handle = Entrez.esearch(db="assembly", term=f'"{species_name}"[Organism] AND "latest refseq"[Filter]')
    assembly_records = Entrez.read(handle)
    handle.close()
    
    if not assembly_records["IdList"]:
        raise ValueError(f"No RefSeq assemblies found for {species_name}")
    
    assembly_id = assembly_records["IdList"][0]
    
    # Get assembly details
    handle = Entrez.esummary(db="assembly", id=assembly_id)
    assembly_summary = Entrez.read(handle)
    handle.close()
    
    # Get the latest chromosomal sequence
    ftp_path = assembly_summary["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"]
    if not ftp_path:
        raise ValueError(f"No RefSeq FTP path found for {species_name}")
    
    # Extract the genomic sequence file name
    basename = ftp_path.split("/")[-1]
    gbff_file = f"{basename}_genomic.gbff.gz"
    gbff_url = f"{ftp_path}/{gbff_file}"
    
    # Download and parse (you might need to use requests or urllib)
    import requests
    import gzip
    from io import BytesIO
    
    response = requests.get(gbff_url)
    if response.status_code == 200:
        with gzip.open(BytesIO(response.content), 'rt') as f:
            genome_record = seq_read(f, "genbank")
        return genome_record
    else:
        raise ValueError(f"Failed to download genome from {gbff_url}")


def fetch_and_save_genome_file(gene_seq, save_dir="genomes"):
    """
    Descarga y guarda el archivo .gbff.gz del genoma RefSeq para la especie dada.
    """
    species_name = gene_seq.annotations["organism"]

    # Buscar en base de datos Assembly
    handle = Entrez.esearch(db="assembly", term=f'"{species_name}"[Organism] AND "latest refseq"[Filter]')
    assembly_records = Entrez.read(handle)
    handle.close()

    if not assembly_records["IdList"]:
        raise ValueError(f"No RefSeq assemblies found for {species_name}")

    assembly_id = assembly_records["IdList"][0]

    # Obtener resumen del ensamblado
    handle = Entrez.esummary(db="assembly", id=assembly_id)
    assembly_summary = Entrez.read(handle)
    handle.close()

    ftp_path = assembly_summary["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"]
    if not ftp_path:
        raise ValueError(f"No RefSeq FTP path found for {species_name}")

    # Crear URL del archivo .gbff.gz
    basename = ftp_path.split("/")[-1]
    gbff_file = f"{basename}_genomic.gbff.gz"
    gbff_url = f"{ftp_path}/{gbff_file}"

    # Crear carpeta si no existe
    os.makedirs(save_dir, exist_ok=True)
    local_file_path = os.path.join(save_dir, gbff_file)

    # Descargar y guardar
    response = requests.get(gbff_url, stream=True)
    if response.status_code == 200:
        with open(local_file_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"Archivo guardado en: {local_file_path}")
        return local_file_path
    else:
        raise ValueError(f"Error descargando {gbff_url} (status code: {response.status_code})")

def describe_genome(gbff_path):
    records = list(SeqIO.parse(gbff_path, "genbank"))
    
    for record in records:
        print(f"ID: {record.id}")
        print(f"Descripción: {record.description}")
        print(f"Longitud de secuencia: {len(record.seq)}")
        print("Anotaciones:")
        for key, value in record.annotations.items():
            print(f"  {key}: {value}")
        print("\nCaracterísticas (features):")
        for feature in record.features:
            print(f" - Tipo: {feature.type}")
            print(f"   Location: {feature.location}")
            if "gene" in feature.qualifiers:
                print(f"   Gene: {feature.qualifiers['gene']}")
            if "product" in feature.qualifiers:
                print(f"   Producto: {feature.qualifiers['product']}")
        print("\n" + "="*50 + "\n")


def fetch_genome_by_sequence(seq_record, save_dir="genomes"):
    """
    A partir de un SeqRecord con anotaciones, obtiene el organismo,
    busca en Assembly la última RefSeq disponible y descarga el genoma completo.
    """
    # 1. Obtener nombre del organismo
    if "organism" not in seq_record.annotations:
        raise ValueError("El SeqRecord no contiene anotación 'organism'")
    organism_name = seq_record.annotations["organism"]
    print(f"Buscando genoma para organismo: {organism_name}")

    # 2. Buscar ensamblados en Assembly para organismo
    search_term = f'"{organism_name}"[Organism] AND "latest refseq"[Filter]'
    handle = Entrez.esearch(db="assembly", term=search_term, retmax=5)
    assembly_search = Entrez.read(handle)
    handle.close()

    if not assembly_search["IdList"]:
        raise ValueError(f"No se encontraron ensamblados RefSeq para {organism_name}")

    assembly_id = assembly_search["IdList"][0]

    # 3. Obtener resumen del ensamblado
    handle = Entrez.esummary(db="assembly", id=assembly_id)
    assembly_summary = Entrez.read(handle)
    handle.close()

    ftp_path = assembly_summary["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"]
    if not ftp_path:
        raise ValueError(f"No se encontró ruta FTP para RefSeq de {organism_name}")

    # 4. Construir URL para descargar archivo .gbff.gz
    basename = ftp_path.split("/")[-1]
    gbff_file = f"{basename}_genomic.gbff.gz"
    gbff_url = f"{ftp_path}/{gbff_file}"

    # 5. Crear directorio para guardar archivo
    os.makedirs(save_dir, exist_ok=True)
    local_path = os.path.join(save_dir, gbff_file)

    # 6. Descargar y guardar archivo
    print(f"Descargando archivo desde {gbff_url} ...")

    if gbff_url.startswith("ftp://"):
        # Usar urllib para FTP
        urllib.request.urlretrieve(gbff_url, local_path)
        print(f"Archivo guardado en {local_path}")
        return local_path
    else:
        # En caso de HTTP(S), usar requests
        import requests
        response = requests.get(gbff_url, stream=True)
        if response.status_code == 200:
            with open(local_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f"Archivo guardado en {local_path}")
            return local_path
        else:
            raise RuntimeError(f"Error descargando el archivo: status code {response.status_code}")