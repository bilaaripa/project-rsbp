import os
import pandas as pd
from flask import Flask, render_template, request
from flask_wtf import FlaskForm
from wtforms import FileField
from werkzeug.utils import secure_filename
import numpy as np
from matplotlib import pyplot as plt
import biotite.sequence.phylo as phylo
import biotite.sequence.graphics as graphics
from datetime import datetime
import plotly.express as px
from dash import Dash, dcc, html, Input, Output, callback
import dash_bio as dashbio
from scipy.cluster.hierarchy import linkage, leaves_list
# import upgma

app = Flask(__name__)
# dash_app = Dash(__name__, server=app, url_base_pathname='/dash/')

dmatrix = []
vlist = []
dendro_leaves = ()
df = 0
visualization_method = "none"
neighbor_joining_tree = 0

# Create View
@app.route('/')
def welcome():
    return render_template("welcome.html")

@app.route('/index')
def index():
    return render_template("index.html")

# nyoba nyoba masih belum bisa
@app.route("/Upload", methods=['POST'])
def upload():
    global visualization_method
    dataset_file = request.files.get('file')
    visualization_method = request.form.get('visualization_method')

    if dataset_file and dataset_file.filename.endswith('.csv'):
        # Dapatkan timestamp untuk disertakan dalam nama file
        timestamp = datetime.now().strftime("%Y%m%d%H%M%S")

        csv_folder = os.path.join("static", "csv")
        os.makedirs(csv_folder, exist_ok=True)  # Buat folder jika belum ada
        csv_path = os.path.join(csv_folder, f'dataset_{timestamp}.csv')
        dataset_file.save(csv_path)

        print(f"File CSV berhasil disimpan di: {csv_path}")
        
        df = pd.read_csv(csv_path)
        varietieslist = df['N'].tolist()
        global vlist 
        vlist = varietieslist
        df.drop('N', axis=1, inplace=True)
        
        # create matrix varieties
        columns   = list(df)
        rows      = list(df.index)

        varieties   =  [];

        for i in rows:
            k = 0
            varieties.append([])

            for j in range(len(columns)):
                if k%2 == 0 :
                    lc1 = df[columns[j]][i]
                    lc2 = df[columns[j+1]][i]
                    varieties[i].append([lc1, lc2])
                k+=1
    
        # create smlist
        ploidy = 2
        locus = 38
        baris = 1
        kolom = 0
        smlist = []
        index = 0

        for i in range(105):
            content = 0
            mllist = []
            mlindex = 0

            if i == 0:
                smlist.append([])

            for j in varieties[kolom]:
                ml = 0
                if j[0] in varieties[baris][content]:
                    ml+=1

                if j[1] in varieties[baris][content]:
                    ml+=1

                mlphy = ml/ploidy
                mllist.append(mlphy)

                content+=1

            totalml = sum(mllist)
            sm = 1-((1/locus)*totalml)
            smlist[index].append(sm)

            baris+=1
            mlindex+=1
            if baris > 14:
                kolom+=1
                baris = kolom+1
                index+=1

                if i < 104:
                    smlist.append([])
        
        # print(smlist)
        # create distancematrix
        distancematrix = []
        global dmatrix
        r = 0
        c = 0
        r1 = 0
        r2 = 0
        c1 = 0
        rnol = 0
        for i in range(15):
            c2 = i
            distancematrix.append([])

            for j in range(15):
                if j == rnol:
                    distancematrix[r].append(0)
                elif j < rnol:
                    if j == 0:
                        distancematrix[r].append(smlist[0][c1])
                        c1+=1
                        c2-=1
                    else:
                        c2-=1
                        distancematrix[r].append(smlist[j][c2])
                else:
                    distancematrix[r].append(smlist[i][c])
                    c+=1

            c=0
            r+=1
            rnol+=1

        dmatrix = distancematrix
        
        # Tentukan path lengkap untuk menyimpan file di /static/img
        if visualization_method == 'UPGMA':
            # Lakukan visualisasi Dendogram UPGMA
            distances = np.array(distancematrix)
            tree = phylo.upgma(distances)  # Menggunakan UPGMA
            fig, ax = plt.subplots(figsize=(7.0, 4.0))
            graphics.plot_dendrogram(ax, tree, labels=varietieslist)
            fig.tight_layout()

            # Tentukan path lengkap untuk menyimpan file di /static/img
            image_path_upgma = f'static/img/dendrogram_upgma_{timestamp}.png'
            fig.savefig(image_path_upgma)

            saved_tree = tree

            plt.clf()  # Clear plot

            print(saved_tree)

            # Visualisasi dash UPGMA
            # Menampilkan Hasil Akhir
            return render_template("vizresult.html", metode_option="UPGMA", image_path=image_path_upgma, tree=saved_tree)

        elif visualization_method == 'neighbor_joining':
            # Lakukan visualisasi Dendogram Neighbor Joining
            distances = np.array(distancematrix)
            tree = phylo.neighbor_joining(distances)
            fig, ax = plt.subplots(figsize=(7.0, 4.0))
            graphics.plot_dendrogram(ax, tree, labels=varietieslist)
            fig.tight_layout()

            # Tentukan path lengkap untuk menyimpan file di /static/img
            image_path_neighbor_joining = f'static/img/dendrogram_neighbor_{timestamp}.png'
            fig.savefig(image_path_neighbor_joining)
            plt.clf()  # Clear plot

            # Visualisasi dash Neighbor Joining
            # Menampilkan Hasil Akhir
            return render_template("vizresult.html", metode_option="Neighbor Joining", image_path=image_path_neighbor_joining)

        else:
            return "Metode visualisasi tidak valid"

    else:
        return "No file uploaded."

app_dash = Dash(__name__, server=app, url_base_pathname='/dashboard/')    

app_dash.layout = html.Div([
    "Rows to display",
    dcc.Dropdown(
        id='my-default-clustergram-input',
        options=[
            {'label': var, 'value': i} for i, var in enumerate(vlist)
        ],
        value=list(range(15)),
        multi=True
    ),
    html.Div(id='my-default-clustergram')
])

@callback(
    Output('my-default-clustergram', 'children'),
    Input('my-default-clustergram-input', 'value')
)

def update_clustergram(rows):
    global dmatrix
    global vlist
    global dendro_leaves
    global visualization_method
    global neighbor_joining_tree
    
    if len(rows) < 2:
        return "Please select at least two rows to display."
    
    if visualization_method == "UPGMA":
        return dcc.Graph(figure=dashbio.Clustergram(
            data=df.loc[vlist].values[dendro_leaves][:, dendro_leaves],  # Apply dendrogram order to the data
            column_labels=[vlist[i] for i in dendro_leaves],  # Reorder column labels
            row_labels=[vlist[i] for i in dendro_leaves],  # Reorder row labels
            color_threshold={
                'row': 250,
                'col': 700
            },
            hidden_labels='row',
            height=800,
            width=700
        ))
        
    elif visualization_method == "neighbor_joining":
        return dcc.Graph(figure=dashbio.Clustergram(
            data=dmatrix,
            column_labels=vlist,
            row_labels=vlist,
            color_threshold={
                'row': 250,
                'col': 700
            },
            hidden_labels='row',
            height=800,
            width=700,
            tree=neighbor_joining_tree
        ))

# Helper function to calculate the neighbor joining tree
def neighbor_joining(dm, names):
    n = len(dm)
    if n != len(names):
        return None

    def min_element(matrix):
        min_val = float("inf")
        min_i, min_j = -1, -1
        for i in range(len(matrix)):
            for j in range(i+1, len(matrix[i])):
                if matrix[i][j] < min_val:
                    min_val = matrix[i][j]
                    min_i, min_j = i, j
        return min_i, min_j

    def add_new_node(dm, node_name, indices):
        n = len(dm)
        new_dm = np.zeros((n + 1, n + 1))
        new_names = names + [node_name]

        for i in range(n):
            for j in range(i + 1, n):
                if i in indices and j in indices:
                    new_dm[n][i] = (n - 2) * dm[i][j] + sum(dm[i][k] for k in indices if k != i) + sum(
                        dm[j][k] for k in indices if k != j)
                    new_dm[i][n] = new_dm[n][i]
                else:
                    new_dm[i][n] = (dm[i][j] + dm[i][j]) / 2
                    new_dm[n][i] = new_dm[i][n]

        return new_dm, new_names

    tree = []
    while n > 2:
        i, j = min_element(dm)
        new_node = f'({names[i]}, {names[j]})'
        tree.append(new_node)

        dm, names = add_new_node(dm, new_node, [i, j])
        n = len(dm)

    tree.append(f'({names[0]}, {names[1]})')
    return tree[0]


@app.route('/dash')
def dash():   
    global vlist
    global dmatrix
    global dendro_leaves
    global df
    global visualization_method
    global neighbor_joining_tree
    
    if visualization_method == "neighbor_joining":
        # Calculate the neighbor joining tree
        neighbor_joining_tree = neighbor_joining(dmatrix, vlist)
        
    elif visualization_method == "UPGMA":
        # Convert the distance matrix to a Pandas DataFrame
        df = pd.DataFrame(dmatrix, columns=vlist, index=vlist)
        
        # Create dendrogram data using UPGMA clustering
        dendro_data = linkage(df.values, method='average')  # UPGMA clustering
        dendro_leaves = leaves_list(dendro_data)
    
    app_dash.layout = html.Div([
        "Rows to display",
        dcc.Dropdown(
            id='my-default-clustergram-input',
            options=[
                {'label': var, 'value': i} for i, var in enumerate(vlist)
            ],
            value=list(range(15)),
            multi=True
        ),
        html.Div(id='my-default-clustergram')
    ])
    
    return app_dash.index()

if __name__ == "__main__":
    app.run(debug=True)
