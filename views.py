from django.shortcuts import render
from .clusters import *
import urllib
from django.core.files.storage import FileSystemStorage
import os
from django.core.cache import cache
from django.http import HttpResponse
import csv
from django.conf import settings
from index.views import color_palette
import random
import string

def index(request):
    return render(request, "clusters/index.html")


def get_structures(request, pdb=None):
    """
    Computes hydrophobic clusters for specific query or pdb id
    """
    #Directly from URL
    if pdb is not None:
        pdb = pdb.strip()
        if len(pdb) != 4:
            print("error on the pdb ")
            return render(request, "clusters/index.html", {
                'error_message': "PDB ID not valid"
            })
        url = f"https://files.rcsb.org/download/{pdb}.pdb"
        try:
            urllib.request.urlopen(url)
        except:
            return render(request, "clusters/index.html", {
                'error_message': "Failed to download protein from the PDB database. "
                                 "Check that the PDB code is correct and try again"
            })
        mol = Molecule(pdb,validateElements=False)

    # A file was uploaded
    elif request.POST and request.FILES:
        myfile = request.FILES['myfile']
        fs = FileSystemStorage()
        filename = fs.save(myfile.name, myfile)
        if filename.endswith('.pdb') or filename.endswith('.ent') or filename.endswith('.pdb.gz'):
            mol = Molecule(os.path.join(settings.MEDIA_ROOT, filename), validateElements=False)
        else:
            return render(request, "clusters/index.html", {
                'error_message': "File extension not supported"
            })

    # A pdb code was specified
    elif request.POST and request.POST['pdbcode']:
        pdb = request.POST['pdbcode']
        pdb = pdb.strip()
        if len(pdb) != 4:
            print("error on the pdb ")
            return render(request, "clusters/index.html", {
                'error_message': "PDB ID not valid"
            })
        url = f"https://files.rcsb.org/download/{pdb}.pdb"
        try:
            urllib.request.urlopen(url)
        except:
            return render(request, "clusters/index.html", {
                'error_message': "Failed to download protein from the PDB database. "
                                 "Check that the PDB code is correct and try again"
            })
        mol = Molecule(pdb,validateElements=False)

    # Something wrong happened because none of the above worked
    else:
        return render(request, "clusters/trial.html",{
                'error_message': "PDB/File not valid"
            }
        )

    clusters = compute_hydrophobic_clusters(mol)
    if not clusters:
        return render(request, "clusters/trial.html",{
                'error_message': "The selected protein does not contain hydrophobic clusters"
            }
        )
    dict_clusters = []
    for cluster in clusters:
        dict_clusters.append({
            "area":cluster.area,
            'residues': cluster.residues,
            'chains': cluster.chains,
            "contacts": cluster.contacts,
            "ratio_contacts_residue": cluster.ratio_contacts_residue,
            "ratio_area_residue": cluster.ratio_area_residue
        })
    error_message = '' if clusters else 'Structure does not contain hydrophobic clusters'

    cache.set("data", dict_clusters, None)

    if pdb:
        return render(request, "clusters/trial.html",{
            'pdb':pdb,
            'clusters': dict_clusters,
            'error_message': error_message,
        })
    else:

        return render(request, "clusters/trial.html",{
            'filename':filename,
            'clusters': dict_clusters,
            'error_message': error_message,
        })

def documentation(request):
    return render(request, "clusters/documentation.html")


def compute_hydrophobic_clusters(mol,
                                 sel: str = "protein and not backbone and noh and resname ILE VAL LEU",
                                 cutoff_area: float = 10):
    """
    :param chain: Chain in the PDB to compute the hydrophobic clusters. Examples: "A", "A B C". Default: "A"
    :param sel: VMD selection on which to compute the clusters. Default is every sidechain heavy atom ILE, VAL and LEU residues. "protein and not backbone and noh and resname ILE VAL LEU"
    :return: A representation for each cluster
    """
    clusters = None

    resids = mol.get("resid", sel="protein and name CA and resname ILE VAL LEU") #ILV cluster residues
    chains = mol.get("chain", sel="protein and name CA and resname ILE VAL LEU") # chains of ILV cluster residues
    dims = len(resids) # length of ILV clusters
    indices = mol.get("index", sel=f"{sel}") #the indices of the atoms in ILV clusters

    # get a dictionary of atom index and resid position in resids list
    atom_to_residposition = {}
    for index in indices:
        resid = mol.get("resid", sel=f"index {index}")[0]
        chain = mol.get("chain", sel=f"index {index}")[0]
        index_residue = [j for j, residue in enumerate(resids) if (residue == resid and chains[j] == chain) ][0]
        atom_to_residposition[index] = index_residue

    contacts = np.zeros((dims, dims))

    for index in indices:
        a = Atom(index, mol)
        if not a.neighbor_indices.any():
            continue
        contacts = fill_matrices(a, mol, contacts, indices, atom_to_residposition)

    graph = create_graph(contacts, resids, chains, cutoff_area=cutoff_area)
    comp, _ = label_components(graph)
    if comp.a.any():
        clusters = add_clusters(graph, comp)
    else:
        print("There are not residues in contact for this selection")

    return clusters

def download_table(request):
    """A view that streams a CSV file."""
    # Generate a sequence of rows. The range is based on the maximum number of
    # rows that can be handled by a single sheet in most spreadsheet
    # applications.


    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="fragment_table.csv"'
    writer = csv.writer(response)

    writer.writerow(['Cluster ID','Area','Number of contacts',
                     'Contacts/Residue','Area/Residue'])
    query_set = cache.get("data")

    output=[]
    for index, cluster in enumerate(query_set):
        output.append([index, round(cluster['area'],2), round(cluster['contacts'],2),
                       round(cluster['ratio_contacts_residue'],2), round(cluster['ratio_area_residue'],2)])

    # CSV Data
    writer.writerows(output)
    return response

def download_session(request, filename=None, pdb=None):
    """
    A view that generates an exportable PyMol session
    """
    print(filename, pdb)
    dict_clusters = cache.get("data")

    #If it comes from a pdb code we need to save the pdb in the media.
    if pdb:
        mol=Molecule(pdb,validateElements=False)
        rndstr= ''.join(random.choice(string.ascii_letters) for i in range(10))
        filename = f'{pdb}_{rndstr}.pdb'
        path_filename = f'{settings.MEDIA_ROOT}/{filename}'
        mol.write(path_filename)

    #define the media url
    full_url = urllib.parse.urljoin('http://proteintools.uni-bayreuth.de', settings.MEDIA_URL)

    # first line of the pymol session
    sessions = f"load {full_url}{filename}, protein\n color white, protein;\n"

    # iterate over clusters and color them
    for index, cluster in enumerate(dict_clusters):
        if len(set(cluster['chains'])) == 1:
            sessions += f"select cluster{index}, chain {cluster['chains'][0]} and resi "
            for resi in cluster['residues']:
                sessions += f"{resi}+"
        else:
            sessions += f"select cluster{index}, "
            for j,resi in enumerate(cluster['residues']):
                if j < len(cluster['residues']) - 1:
                    sessions += f"(chain {cluster['chains'][j]} and resi {resi}) or "
                else: # last pair chain-residue, remove last OR.
                    sessions += f"(chain {cluster['chains'][j]} and resi {resi})"
        sessions += f"\n"
        sessions += f"show spheres, cluster{index}\n"
        sessions += f"set_color {color_palette[index][0]}, {color_palette[index][1:4]}  \n"
        sessions += f"color {color_palette[index][0]}, cluster{index}\n"

    sessions += f"set cartoon_color, white\n"

    response = HttpResponse(sessions, content_type='application/x-pymol')

    response['Content-Disposition'] = 'attachment; filename=session.pml'

    return response