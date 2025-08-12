import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import umap
import matplotlib.pyplot as plt
from sklearn.utils import shuffle
from tqdm import tqdm
from matplotlib.colors import ListedColormap
from scipy.spatial import ConvexHull
from matplotlib.patches import Polygon
from multiprocessing import Pool, cpu_count
import joblib
from sklearn.manifold import TSNE
import seaborn as sns
 
# --- CONFIG ---
CSV_FILES = [
    "/home/qg622/qianrong/Y2/ScaffoldSplitsOverestimateVS/data/clustering_id_k7.csv",
    "/home/qg622/qianrong/Y2/ScaffoldSplitsOverestimateVS/figures/real-world_diversity/fda_ingredients_smiles.csv",
    "/home/qg622/qianrong/Y2/ScaffoldSplitsOverestimateVS/figures/real-world_diversity/zinc.csv"
    ]
CHUNKSIZE = 100000
MAX_CHUNKS = 300
N_BITS = 256
RADIUS = 2
DATASET_NAMES = ["NCI60", "FDA Approved", "ZINC"]
N_PROCESSES = 100
 
smiles_list = []
embedding_list = []
labels = []
 
# --- Helper ---
def mol_to_fp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=RADIUS, nBits=N_BITS)
    return np.array(fp)
 
def process_smiles_chunk(args):
    smiles_chunk, label = args
    result_fps = []
    result_labels = []
    for smi in smiles_chunk:
        fp = mol_to_fp(smi)
        if fp is not None:
            result_fps.append(fp)
            result_labels.append(label)
    return result_fps, result_labels
 
# --- UMAP reducer ---
use_tsne = True  # Set to True to use t-SNE instead of UMAP
if use_tsne:
    reducer = TSNE(n_components=2, random_state=42, perplexity=30, n_jobs=-1)
    method_label = 't-SNE'
else:
    reducer = umap.UMAP(n_components=2, random_state=42, verbose=True)
    method_label = 'UMAP'
# reducer = umap.UMAP(n_components=2, random_state=42)
reducer_fitted = False
 
# --- Read, process, and reduce in chunks with multiprocessing ---
if method_label == 't-SNE':
    all_fps, all_lbls = [], []
 
with Pool(processes=N_PROCESSES) as pool:
    for idx, file in enumerate(CSV_FILES):
        total_count = 0
        chunks = pd.read_csv(file, chunksize=CHUNKSIZE)
        for i, chunk in enumerate(tqdm(chunks, desc=f"Processing {file}", unit="chunk")):
            if i >= MAX_CHUNKS:
                break
            smiles_col_name = 'smiles' if 'smiles' in chunk.columns else 'SMILES'
            smiles_col = chunk[smiles_col_name].dropna().tolist()
            split_chunks = np.array_split(smiles_col, N_PROCESSES)
            results = pool.map(process_smiles_chunk, [(sublist, idx) for sublist in split_chunks])
            fps, lbls = [], []
            for fps_chunk, lbl_chunk in results:
                fps.extend(fps_chunk)
                lbls.extend(lbl_chunk)
            total_count += len(fps)
            labels.extend(lbls)
            if not fps:
                continue
            fps_array = np.array(fps)
            if method_label == 't-SNE':
                all_fps.extend(fps_array.tolist())
                all_lbls.extend(lbls)
            else:
                if not reducer_fitted:
                    embedding_chunk = reducer.fit_transform(fps_array)
                    reducer_fitted = True
                else:
                    embedding_chunk = reducer.transform(fps_array)
                embedding_list.append(embedding_chunk)
        print(f"{file}: {total_count} molecules processed.")
 
if method_label == 't-SNE' and all_fps:
    all_fps_array = np.array(all_fps)
    embedding_list = reducer.fit_transform(all_fps_array)
    labels = all_lbls
    
# --- Final Embedding ---
embedding = np.vstack(embedding_list)
labels = np.array(labels)
embedding, labels = shuffle(embedding, labels, random_state=42)
 
# --- Save embedding and labels ---
joblib.dump((embedding, labels), f"embedding_labels_{method_label.lower()}.pkl")
 
# --- Plot ---
reds = plt.cm.get_cmap('Reds', 7)
cmap = ListedColormap(reds(np.linspace(0.3, 1.0, 3)))
fig, ax = plt.subplots(figsize=(10, 10))
sns.set_style('whitegrid')
 
for i in range(2, -1, -1):
    idx = labels == i
    color = cmap(i)
    ax.scatter(embedding[idx, 0], embedding[idx, 1], s=5, label=DATASET_NAMES[i], alpha=0.5, color=color)
    points = embedding[idx]
    if len(points) >= 3:
        hull = ConvexHull(points)
        polygon = Polygon(points[hull.vertices], edgecolor=color, fill=False, linewidth=1)
        ax.add_patch(polygon)
 
# ax.set_title(f"{method_label} of Morgan Fingerprints from 3 SMILES Datasets")
# ax.set_xlabel(f"{method_label}-1")
# ax.set_ylabel(f"{method_label}-2")
ax.legend(fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
plt.tight_layout()
plt.savefig(f'{method_label.lower()}_real_world_distribution.png',dpi=600)
plt.show()
 
 