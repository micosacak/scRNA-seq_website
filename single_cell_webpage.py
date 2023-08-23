import os
import string, difflib
import _pickle as cPickle
from flask import Flask, render_template, request, redirect, url_for
app = Flask(__name__)
app.config['SECRET_KEY'] = "????????????????????????????????"
global tsne_image_files
global intrxs_image_files 
global ensemblid2genename
def getPossGenes(my_org, gene_name, max_len = 3, i = 1):
	# a naive algorithm to search for the gene names.
	indexes = getIndexes(my_org,gene_name[0:i])
	while len(getIndexes(my_org,gene_name[0:i])) >= max_len and i <= len(gene_name):
		indexes = getIndexes(my_org,gene_name[0:i])
		i += 1
	possible_genes = [my_org[idx] for idx in indexes]  
	possible_genes.sort()
	return possible_genes[0:50] 
def checkword(word):
	word = word.replace(" ","")
	while word.startswith("."):
		word = word[1:]
	while word.endswith("."):
		word = word[:-1]
	word = word.replace(":","_")
	return word
def getIndexes(my_list, my_char):
    return [i for i in range(0, len(my_list)) if my_list[i].startswith(my_char) or my_list[i].startswith(my_char.lower())]
def get_image_names(target_file = None):
    my_folders = os.listdir(f"./static/image_files/{target_file}/")
    file_names = {}
    print(f"./static/image_files/{target_file}/")
    for my_folder in my_folders:
        file_names[my_folder] = []
        file_names[my_folder] = os.listdir(f"./static/image_files/{target_file}/{my_folder}/")
        file_names[my_folder] = [file_name[:-4] for file_name in file_names[my_folder] if file_name.endswith(".pdf")]
    return file_names
def get_end2gene():
	input = open("static/genes.tsv")
	map = {}
	for line in input:
		ensid, geneid = line.strip("\n").split("\t")
		map[ensid] = geneid
	return map
if "data" not in os.listdir("static"):
	os.mkdir("static/data")
if("tsne_file_names.p" in os.listdir("./static/data/")):
	tsne_file_names = cPickle.load(open("./static/data/tsne_file_names.p","rb"))
else:
	tsne_file_names = get_image_names(target_file = "TSNE")
	cPickle.dump(tsne_file_names, open("./static/data/tsne_file_names.p","wb"))
if("intrxs_file_names.p" in os.listdir("./static/data/")):
	intrxs_file_names = cPickle.load(open("./static/data/intrxs_file_names.p","rb"))
else:
	intrxs_file_names = get_image_names(target_file = "InteractionsByGene")
	cPickle.dump(intrxs_file_names, open("./static/data/intrxs_file_names.p","wb"))
if("ens2gene.p" in os.listdir("./static/data/")):
	ens2gene = cPickle.load(open("./static/data/ens2gene.p","rb"))
else:
	ens2gene = ens2gene = get_end2gene()
	cPickle.dump(ens2gene, open("./static/data/ens2gene.p","wb"))
tsne_image_files = tsne_file_names
intrxs_image_files = intrxs_file_names
ensemblid2genename = ens2gene
@app.route("/")
@app.route("/home/")
def home():
    return render_template("home.html", title="Home")
@app.route("/DisclaimerPage/", methods = ['GET','POST'])
def DisclaimerPage():
	return render_template("DisclaimerPage.html", title="Disclaimer")
@app.route("/TSNE/")
def get_TSNE_plots():
    return render_template("TSNE.html",title="TSNE Plots", first_selections = [{"name":"ALL Cells"},{"name":"Progenitor Cells"},{"name":"Neurons"}])
@app.route("/TSNE/showTSNE/", methods = ['GET','POST'])
def show_TSNE_plots():    
    sample_ids = {"ALL Cells":"ALL", "Progenitor Cells":"PC","Neurons":"Neuron"}
    sample_id = request.form.get('sample_id_select', )
   # gene_name = checkword(request.form["gene_name_select"])  # updated on 13.08.2019
    gene_name = request.form["gene_name_select"] # the above line will cause bug if the gene_name has ":" in. # updated on 13.08.2019
    if gene_name.startswith("ENSDARG"):
        if gene_name in ensemblid2genename:
            gene_name = ensemblid2genename[gene_name]
    # check the gene_name to remove space " ", ".", or ":". #updated on 13.08.2019
    gene_name = checkword(gene_name)  # added on 13.08.2019
    if (gene_name not in tsne_image_files[f"{sample_ids[sample_id]}_Ftr"]) and (gene_name.lower() not in tsne_image_files[f"{sample_ids[sample_id]}_Ftr"]):
        if len(gene_name) <= 1:
            return render_template("showTSNEerror.html", title="TSNE Plots",msgs = f"The gene name is too short!", possible_genes = ["No Suggestion!"], sample_id = sample_ids[sample_id])
        else:
            # a naive algorithm to search for the gene names.
            my_org = tsne_image_files[f"{sample_ids[sample_id]}_Ftr"]
            possible_genes = getPossGenes(my_org, gene_name)			
            return render_template("showTSNEerror.html",title="TSNE Plots", msgs = f"Could not find the gene: {gene_name}", possible_genes = possible_genes, sample_id = sample_ids[sample_id])
    else:
        if gene_name.lower() in tsne_image_files[f"{sample_ids[sample_id]}_Ftr"]:
            gene_name = gene_name.lower()
        return render_template("showTSNE.html", title="TSNE Plots",sample_id = sample_ids[sample_id], gene_name = gene_name)
@app.route("/GO/")
def get_GO_plots():
    return render_template("GO.html",title="GO Plots", first_selections = [{"name":"AB42_vs_PBS"},{"name":"IL4_vs_PBS"}], second_selections = [{"name":"Biological Process"},{"name":"Molecular Function"},{"name":"Cellular Component"}, {"name":"KEGG Pathway"}],third_selections = [{"name":"Im"},{"name":"NN_0"},{"name":"NN_1"}, {"name":"NN_2"},{"name":"NN_3"},{"name":"NN_4"},{"name":"NN_5"},{"name":"NN_6"},{"name":"NN_7"},{"name":"OPCOD_1"},{"name":"OPCOD_2"},{"name":"PC_0"},{"name":"PC_1"},{"name":"PC_2"},{"name":"PC_3"},{"name":"PC_4"},{"name":"PC_5"},{"name":"PC_6"},{"name":"PC_7"},{"name":"PC_8"}])
@app.route("/GO/showGO/", methods = ['GET','POST'])
def show_GO_plots():    
    gopath_ids = {"Biological Process":"BP_Up", "Molecular Function":"MF_Up","Cellular Component":"CC_Up","KEGG Pathway":"kegg_Up"}
    comparison_id = request.form.get('comparison_id_select')
    gopath_id = request.form.get('gopath_id_select')
    cellid = request.form["cell_id_select"]
    return render_template("showGO.html", title="GO Plots",gopath_id = gopath_ids[gopath_id],comparison_id = comparison_id, cellid = cellid)
@app.route("/InteractionMaps/")
def InteractionMaps():
	return render_template("InteractionMaps.html", title="Interaction Maps")
@app.route("/InteractionMaps/ByPath/")
def get_InxMapByPath():
    return render_template("ByPath.html",title="Interaction Maps", first_selections = [{"name":"AB42_Only"},{"name":"IL4_Only"},{"name":"PBS_Only"},{"name":"AB42_vs_PBS"},{"name":"IL4_vs_PBS"}], second_selections = [{"name":"agrn"},{"name":"appa"},{"name":"ctgfa"},{"name":"cxc"},{"name":"edil3a"},{"name":"fgf"},{"name":"gnai"},{"name":"hbegf"},{"name":"notch"},{"name":"ptn"},{"name":"serpine1"},{"name":"serpine2"}])
@app.route("/InteractionMaps/ByPath/showInxMapByPath/", methods = ['GET','POST'])
def show_InxMapByPath():    
    comparison_id = request.form.get('comparison_id_select')
    pathname_id = request.form.get('pathname_id_select')
    return render_template("showInxMapByPath.html", title="Interaction Maps",pathname_id = pathname_id,comparison_id = comparison_id)   
@app.route("/InteractionMaps/ByGene/")
def get_InxMapByGene():
    return render_template("ByGene.html",title="Interaction Maps", first_selections = [{"name":"AB42_Only"},{"name":"IL4_Only"},{"name":"PBS_Only"},{"name":"AB42_vs_PBS"},{"name":"IL4_vs_PBS"}])
@app.route("/InteractionMaps/ByGene/showInxMapByGene/", methods = ['GET','POST'])
def show_InxMapByGene(): 
    comparison_ids = {"AB42_Only":"AB42_Only","IL4_Only":"IL4_Only","PBS_Only":"PBS_Only","AB42_vs_PBS":"AB42_vs_PBS","IL4_vs_PBS":"IL4_vs_PBS"}
    comparison_id = request.form.get('comparison_id_select')
    # gene_name = checkword(request.form["gene_name_select"]) #updated on 13.08.2019
    gene_name = request.form["gene_name_select"] # the above line will cause  bug if the gene name has ":" in. # updated on 13.08.2019
    if gene_name.startswith("ENSDARG"):
        if gene_name in ensemblid2genename:
            gene_name = ensemblid2genename[gene_name]
    # check the gene_name to remove space " ", "." or ":".
    gene_name = checkword(gene_name) # added on 13.08.2019
    if (gene_name not in intrxs_image_files[f"{comparison_ids[comparison_id]}"]) and (gene_name.lower() not in intrxs_image_files[f"{comparison_ids[comparison_id]}"]):
        if len(gene_name) <= 1:
            return render_template("showInterxsError.html",title="Interaction Maps", msgs = f"The gene name is too short!", possible_genes = ["No Suggestion!"], comparison_id = comparison_ids[comparison_id])
        else:
            # a naive algorithm to search for the gene names.
            my_org = intrxs_image_files[f"{comparison_ids[comparison_id]}"]
            possible_genes = getPossGenes(my_org, gene_name)
            return render_template("showInterxsError.html", title="Interaction Maps", msgs = f"Could not find the gene: {gene_name}", possible_genes = possible_genes, comparison_id = comparison_ids[comparison_id])
    else:
        if gene_name.lower() in intrxs_image_files[f"{comparison_ids[comparison_id]}"]:
            gene_name = gene_name.lower()
            return render_template("showInxMapByGene.html", title="Interaction Maps",gene_name = gene_name,comparison_id = comparison_id) 
@app.route("/help/")
def help():
    return render_template("help.html", title="Help")
@app.route("/Impressum/")
def Impressum():
	return render_template("Impressum.html", title="Impressum")
@app.route("/NewsUpdates/")
def NewsUpdates():
	return render_template("NewsUpdates.html",title="updates")
if __name__ == "__main__":
    app.run(debug=True)
