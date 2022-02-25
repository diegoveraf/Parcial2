#LLamado a las librerias que se van a utilizar para la resolucion del examen
import re
import csv
import pandas as pd
import itertools
from Bio import Entrez


#Funcion de busqeda con el keyword
def download_pubmed(keyword): 
    """En esta función se llama a Entrez para hacer una busqueda en NCBI, conectandose al servidos sin la necesidad de descarar los datos."""
    #"Ecuador genomics [Title/Abstract]"
    Entrez.email = "diego.vera@est.ikiam.edu.ec"
    handle = Entrez.esearch (db = "pubmed",
                             term = keyword,
                             usehistory = "y")
    record = Entrez.read (handle)
    id_list = record ["IdList"]
    webenv = record["WebEnv"]
    query_key = record["QueryKey"]
    handle = Entrez.efetch(db="pubmed", 
                           rettype="medline", 
                           retmode="text", 
                           retstart=0, 
                           retmax = 1500, 
                           webenv = webenv, 
                           query_key = query_key)
    #Lectura de los archivos de nube, y se cambia los saltos de linea
    key1 = handle.read()
    key2 = re.sub(r'\n\s{6}', ' ', key1)
    return (key2)

def mining_pubs(tipo,key2):
    """La función mining_pubs, posee 2 argumentos que corresponden a la lectura de nube y el tipo que es: DP, AU o AD.
    Si el tipo es "DP" recupera el año de publicación del artículo. El retorno es un dataframe con el PMID y el DP_year.
    Si el tipo es "AU" recupera el número de autores por PMID. El retorno es un dataframe con el PMID y el num_auth.
    Si el tipo es "AD" recupera el conteo de autores por país. El retorno es un dataframe con el country y el num_auth."""
    
    #Si el tipo es "PMDI" recupera el año de publicación del artículo. El retorno es un dataframe con el PMID y el DP_year
    PMID = re.findall("PMID- (\d*)", key2) 
    #Si el tipo es "DP" recupera el año de publicación del artículo. El retorno es un dataframe con el PMID y el DP_year
    year = re.findall("DP\s{2}-\s(\d{4})", key2) ##extrae el año de publicación de los articulos 
    #Resultado de la busqueda y dataframe
    numeroYearsPMID = pd.DataFrame()
    numeroYearsPMID["PMID"] = PMID
    numeroYearsPMID["Years"] = year
    
    #Si el tipo es "AU" recupera el número de autores por PMID. El retorno es un dataframe con el PMID y el num_auth. 
    #Patron de busqueda para los autores 
    autoresPapers = key2.split("PMID- ")
    autoresPapers.pop(0) #Se purgan datos
    #Lista para el numero de autores
    numeroAutores = []
    #Loop para contar el numero de autores
    for i in range(len(autoresPapers)):
        num = re.findall("AU -", autoresPapers[i])
        n = (len(num))
        numeroAutores.append(n)
    #Resultado de la busqueda y dataframe
    numeroAutoresPMID = pd.DataFrame()
    numeroAutoresPMID["PMID"] = PMID 
    numeroAutoresPMID["Numero de autores"] = numeroAutores
    
    #Si el tipo es "AD" recupera el conteo de autores por país. El retorno es un dataframe con el country y el num_auth
    #Listas
    paisesAD = []
    #Se busca la informacion esperada
    for line in key2.splitlines():
        if line.startswith("AD  -"):
            paisesAD.append(line[:])
    #Listas
    pubmedPaises1 = []
    pubmedPaises2 = []
    pubmedPaises3 = []
    pubmedPaises4 = []
    pubmedPaises5 = []
    pubmedPaises6 = []
    pubmedPaises7 = []
    pubmedPaisesTotal = []
    for line in key2.splitlines():
        if line.startswith("AD  -"):
            #Se unen todos los ADs en una sola lista
            paisesAD = line[:]
            #Patrones mas repetidos para paises de la busqueda
            #Se buscaron paises que tiene antews numeros y letras o que finalizan en correos
            pubmedPais1 = re.findall(r'\,\s(\w+)\.', paisesAD)
            pubmedPais2 = re.findall(r'\,\s(\w+[^0-9\,]\s\w+[^0-9])\.', paisesAD)
            pubmedPais3 = re.findall(r'\,\s(\w+)\.\s[a-z0-9_\.-]+@[\da-z\.-]+\.[a-z\.]{2,6}', paisesAD)
            pubmedPais4 = re.findall(r'\,\s(\w+[^0-9\,]\s\w+[^0-9])\.\s[a-z0-9_\.-]+@[\da-z\.-]+\.[a-z\.]{2,6}', paisesAD)
            pubmedPais5 = re.findall(r'\,\s(\w+)\. Electronic address:\s[a-z0-9_\.-]+@[\da-z\.-]+\.[a-z\.]{2,6}', paisesAD)
            pubmedPais6 = re.findall(r'\,\s(\w+[^0-9\,]\s\w+[^0-9])\. Electronic address:\s[a-z0-9_\.-]+@[\da-z\.-]+\.[a-z\.]{2,6}', paisesAD)
            pubmedPais7 = re.findall(r'\,\s\w+[0-9\-]\,\s(\w+)\.\n', paisesAD)
            pubmedPaises1.append(pubmedPais1)
            pubmedPaises2.append(pubmedPais2)
            pubmedPaises3.append(pubmedPais3)
            pubmedPaises4.append(pubmedPais4)
            pubmedPaises5.append(pubmedPais5)
            pubmedPaises6.append(pubmedPais6)
            pubmedPaises7.append(pubmedPais7)
    #Se unene las listas de los patrones
    pubmedPaisesTotal = pubmedPaises1 + pubmedPaises2 + pubmedPaises3 + pubmedPaises4 + pubmedPaises5 + pubmedPaises6 + pubmedPaises7
    #A la lista se le saca de los parenteisis y se la ubica en una sola para que el dataframe se pueda usar
    #Luego se hace lo mismo que la tarea de Map of Science
    pubmedPaisesTotal = list(itertools.chain.from_iterable(pubmedPaisesTotal))
    len(pubmedPaisesTotal)
    unicopubmedPaisesTotal = list(set(pubmedPaisesTotal))
    unicopubmedPaisesTotal.sort()
    len(unicopubmedPaisesTotal)
    #Se asigana los datos al diccionario provenientes de las coordendas del mundo
    #Declarar diccionario vacio y abre el archivo de coordenadas
    coordenadasDelMundo = {}
    with open('coordenadas-del-mundo.txt') as f:
        csvr = csv.DictReader(f)
        for row in csvr:
            coordenadasDelMundo[row['Name']] = [row['Latitude'], row['Longitude']]
        #Se crea una lista por cada elemento de unico y se compara con el documento de coordenadas para asi graficarlo y contarlo
    paisesTotalEnOrden = []
    longitudDePaisesDelMundo = []
    latitudDePaisesDelMundo = []
    contDePaisesDelMundo = []
    for i in unicopubmedPaisesTotal:
        if i in coordenadasDelMundo.keys():
            paisesTotalEnOrden.append(i)
            latitudDePaisesDelMundo.append(float(coordenadasDelMundo[i][0]))
            longitudDePaisesDelMundo.append(float(coordenadasDelMundo[i][1]))
            contDePaisesDelMundo.append(pubmedPaisesTotal.count(i))
    #Dataframe de paises
    numeroCountryAU = pd.DataFrame()
    numeroCountryAU["country"] = paisesTotalEnOrden 
    numeroCountryAU["num_auth"] = contDePaisesDelMundo
    
    #Dato que retorna del tipo
    if tipo == 'AD':
        return numeroYearsPMID
    if tipo == 'AU':
        return numeroAutoresPMID
    if tipo == 'PD':
        return numeroCountryAU