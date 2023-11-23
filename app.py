import os
import pandas as pd
from flask import Flask, render_template, request
from flask_wtf import FlaskForm
from wtforms import FileField
from werkzeug.utils import secure_filename
from varietas import varietas
import numpy as np
from matplotlib import pyplot as plt
import biotite.sequence.phylo as phylo
import biotite.sequence.graphics as graphics
from datetime import datetime

app = Flask(__name__)

# Create View
@app.route('/')
def index():
    return render_template("index.html")

# nyoba nyoba masih belum bisa
@app.route("/Upload", methods=['POST'])
def upload():
    dataset_file = request.files.get('file')
    visualization_method = request.form.get('visualization_method')

    print(visualization_method)

    if dataset_file and dataset_file.filename.endswith('.csv'):
        # Dapatkan timestamp untuk disertakan dalam nama file
        timestamp = datetime.now().strftime("%Y%m%d%H%M%S")

        filename = secure_filename(dataset_file.filename)
        dataset_file.save(filename)

        df = pd.read_csv(filename)
        varietieslist = df['N'].tolist()
        df.drop('N', axis=1, inplace=True)
        
        os.remove(filename)
        
        # create matrix varieties
        columns = list(df)
        rows = list(df.index)
        varieties = [[] for _ in rows]

        for i in rows:
            k = 0

        for j in range(len(columns)):
            if k % 2 == 0:
                lc1 = df[columns[j]][i]
                lc2 = df[columns[j + 1]][i]
                varieties[i].append([lc1, lc2])
            k += 1
    
        # create smlist
        ploidy = 2
        locus = 38
        baris = 1
        kolom = 0
        smlist = []
        index = 0
        matrixrange = 2 * locus
        
        for i in range(105):
            content = 0
            mllist = []
            mlindex = 0

            if i == 0:
                smlist.append([])

            for j in varieties[kolom]:
                ml = 0
                if j[0] in varieties[baris][content]:
                    ml += 1

                if j[1] in varieties[baris][content]:
                    ml += 1

                mlphy = ml / ploidy
                mllist.append(mlphy)

                content += 1

            totalml = sum(mllist)
            sm = 1 - ((1 / locus) * totalml)
            smlist[index].append(sm)

            baris += 1
            mlindex += 1
            if baris > 14:
                kolom += 1
                baris = kolom + 1
                index += 1

                if i < 104:
                    smlist.append([])
        
        # create distancematrix
        distancematrix = []
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
                        c1 += 1
                        c2 -= 1
                    else:
                        c2 -= 1
                        distancematrix[r].append(smlist[j][c2])
                else:
                    distancematrix[r].append(smlist[i][c])

            c = 0
            r += 1
            rnol += 1

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
            plt.clf()  # Clear plot

            # Visualisasi dash UPGMA
            # Menampilkan Hasil Akhir
            return render_template("vizresult.html", metode_option="UPGMA", image_path=image_path_upgma)

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

if __name__ == "__main__":
    app.run(debug=True)
