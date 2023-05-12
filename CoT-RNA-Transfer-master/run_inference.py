from model.model import CoT_RNA_Transfer
import argparse
from create_dataset import *
from misc import *

logger = logging.getLogger(__name__)

np.set_printoptions(threshold=sys.maxsize)
def msa_to_embed(msa_path, max_seqs=500, AminoAcids='HETL'):
    tmp_path = msa_path.replace('.fasta', '.fa')
    lines = []
    for line in open(msa_path):
        line = line.strip()
        if not line.startswith(">"):
            new_line = ''
            for l in line:
                if l == 'A':
                    new_line += AminoAcids[0]
                elif l == 'U':
                    new_line += AminoAcids[1]
                elif l == 'C':
                    new_line += AminoAcids[2]
                elif l == 'G':
                    new_line += AminoAcids[3]
                else:
                    new_line += '-'
            lines.append(new_line)
        else:
            lines.append(line)

    if max_seqs is not None:
        lines = lines[:2*max_seqs]     ### 2x for name and sequence

    with open(tmp_path, 'w') as f:
        for line in lines:
            f.write(f"{line}\n")

    if lines[0].startswith(">"):
        L = len(lines[1].strip())

    program = [
        os.path.join(os.path.dirname(__file__), "bin/a3m_to_feat"),
        "--input",
        tmp_path,
        "--max_gap",
        "7",
        "--max_keep",
        "50000",
        "--sample_ratio",
        "1.0",
    ]
    process = subprocess.run(program, capture_output=True)
    assert process.returncode == 0, "Invalid A3M file"
    x = np.copy(np.frombuffer(process.stdout, dtype=np.int8))
    x = x.reshape((-1, L, 7 * 2 + 3)).transpose((0, 2, 1))
    assert (x < 23).all(), "Internal error"
    seq = x[0][0]

    os.remove(tmp_path)

    # return {
    #     "seq": torch.tensor(seq).long()[None].cuda(),
    #     "msa": torch.tensor(x).long()[None].cuda(),
    #     "index": torch.arange(seq.shape[0]).long()[None].cuda(),
    # }

    return {
        "seq": torch.tensor(seq).long()[None],
        "msa": torch.tensor(x).long()[None],
        "index": torch.arange(seq.shape[0]).long()[None],
    }


def main():

    parser = argparse.ArgumentParser()
    # parser.add_argument('--input_MSA', default='RNA_TESTSET/MSA_pydca/RF00500.faclean', type=str)
    parser.add_argument('--input_MSA', default='/home/ubuntu/deepbreaks/data/sars_top500.fasta', type=str)
    parser.add_argument('--model', default='pretrained_models/model.chk', type=str)
    args = parser.parse_args()

    ### print args
    hparams_dict = dict()
    for arg in vars(args):
        hparams_dict[arg] = getattr(args, arg)
        print(arg, getattr(args, arg))

    ### model definition
    model = CoT_RNA_Transfer()

    ### load params
    weight_path = os.path.join(os.path.dirname(__file__), args.model)
    state_dict = torch.load(weight_path)
    model.load_state_dict(state_dict)

    ### move model to GPU
    # model = model.cuda()

    ### traslate MSA from nucletide to amino acids
    adapted_msa = msa_to_embed(args.input_MSA)

    ### evaluate model
    model.eval()
    with torch.no_grad():
        pred, feat = model(adapted_msa)
    pred = pred.cpu()

    L = pred.shape[0]
    mask = torch.full((L, L), -10000)
    for i in range(L):
        for j in range(L):
            if abs(i - j) > 4:
                mask[i, j] = 0
            else:
                pass

    pred = pred.cpu() + mask
    delta = torch.randn(L, L) * 1e-7
    pred = pred + delta + delta.T

    ### save raw output
    dist = pred
    np.savetxt('outputs/dist.txt', dist.numpy())

    ### save top-L prediction
    topk_values, _ = pred.reshape(-1).topk(k=int(2 * 1 * L))
    topk_value = topk_values[-1]
    pred[pred < topk_value] = -10000
    pred[pred >= topk_value] = 1
    pred[pred <= 0] = 0

    np.savetxt('outputs/pred.txt', pred.numpy().astype(int), fmt='%d', delimiter=",")


    # Add.......................................................................03/06/2023

    pre = pred.numpy().astype(int)
    ### get ground truth for input msa

    test_pdb_data_pickle_file = 'RNA_TESTSET_PDB_DATA.pickle'
    if os.path.exists(test_pdb_data_pickle_file):
        with open(test_pdb_data_pickle_file, 'rb') as handle:
            test_pdb_data = pickle.load(handle)
    rna_fam_name = args.input_MSA.split('/')[-1].split('.')[0]
    if rna_fam_name in test_pdb_data:
        test_label = np.ones((L, L)) * -100  ##### ignore index is -100
        for k, v in test_pdb_data[rna_fam_name].items():
            i, j = k[0], k[1]
            if abs(i - j) > 4:
                lbl = distance_to_37(v[-1])
                if lbl <= -1:
                    lbl2 = 100
                elif lbl < 16:
                    lbl2 = 1  ##### lbl is 1 (contact) if distance is smaller than 10A, which corresponds to label 0,1,2,...,15
                elif lbl >= 16:
                    lbl2 = 0  ##### lbl is 0 (non-contact) if distance is larger than 10A, which corresponds to label 16,17,18,....
                test_label[i, j] = lbl2
                test_label[j, i] = lbl2

        test_label[test_label == -100] = 0

    b = test_label
    import matplotlib.pyplot as plt

    def counts(M):

        adjacency_matrix = M

        # Initialize the number of cuts to zero
        num_cuts = 0

        # Iterate over all pairs of vertices in the graph
        for i in range(adjacency_matrix.shape[0]):
            for j in range(i + 1, adjacency_matrix.shape[0]):

                # Check if there is an edge between vertices i and j
                if adjacency_matrix[i, j] == 1:

                    # Iterate over all other pairs of vertices in the graph
                    for k in range(adjacency_matrix.shape[0]):
                        for l in range(k + 1, adjacency_matrix.shape[0]):

                            # Check if there is an edge between vertices k and l
                            if adjacency_matrix[k, l] == 1:

                                # Check if there is a cut between vertices i and j and vertices k and l
                                if (i <= k <= j < l) or (k < i <= l <= j):
                                    num_cuts += 1

        print("Number of cuts:", (num_cuts))

    def diff(n, t):
        z = n
        # t = t.numpy()
        for i in range(len(n)):
            for j in range(len(n)):
                if n[i, j] == 1:
                    if t[i, j] == 0:
                        z[i, j] = 100     #wrong prediction
        z = np.where(z == 1, 50, z)
        cmap = plt.cm.get_cmap('RdYlBu', 3)
        plt.imshow(z, cmap=cmap)
        plt.xticks(np.arange(z.shape[1]))
        plt.yticks(np.arange(z.shape[0]))
        plt.show()
        return z

    z = diff(pre,b)

    # numpy_matrix = np.array([[0, 0, 0, 1], [0, 0, 1, 1], [0, 1, 0, 1], [1, 1, 1, 0]])
    # teb_matrix = np.array([[1, 0, 0, 1], [0, 1, 1, 1], [0, 1, 1, 0], [1, 1, 0, 1]])

    # ...........................................................................
    z1 = np.where(z == 50, 1, z)      # switch 50 back to one
    adjacency_matrix = z1

    # Initialize empty lists to record the indices of the 1s and 100s
    ones_indices = []
    hundreds_indices = []

    # Iterate over the adjacency matrix and record the indices of the 1s and 100s
    for i in range(adjacency_matrix.shape[0]):
        for j in range(adjacency_matrix.shape[1]):
            if adjacency_matrix[i, j] == 1:
                ones_indices.append((i, j))
            elif adjacency_matrix[i, j] == 100:
                hundreds_indices.append((i, j))

    # Initialize a dictionary to record the cuts for each point
    cut_ones = {}
    cut_huns = {}
    # Iterate over the ones and hundreds indices and check for cuts
    for index in ones_indices:
        i, j = index

        # Check if there is an edge between vertices i and j
        if adjacency_matrix[i, j] == 1 or adjacency_matrix[i, j] == 100:

            # Iterate over all other pairs of vertices in the graph
            for k in range(adjacency_matrix.shape[0]):
                for l in range(k + 1, adjacency_matrix.shape[0]):

                    # Check if there is an edge between vertices k and l
                    if adjacency_matrix[k, l] == 1 or adjacency_matrix[k, l] == 100:

                        # Check if there is a cut between vertices i and j and vertices k and l
                        if (i <= k <= j < l) or (k < i <= l <= j):

                            # Record the cuts for each point
                            if index not in cut_ones:
                                cut_ones[index] = 1
                            else:
                                cut_ones[index] += 1

    # Print the cuts for each point
    for point, num_cuts in cut_ones.items():
        print(f"Point {point}: {num_cuts} cuts")
    cutones_array = np.array(list(cut_ones.values()))
    avg_correct = sum(cutones_array) / len(cutones_array)

    for index in hundreds_indices:
        i, j = index

        # Check if there is an edge between vertices i and j
        if adjacency_matrix[i, j] == 1 or adjacency_matrix[i, j] == 100:

            # Iterate over all other pairs of vertices in the graph
            for k in range(adjacency_matrix.shape[0]):
                for l in range(k + 1, adjacency_matrix.shape[0]):

                    # Check if there is an edge between vertices k and l
                    if adjacency_matrix[k, l] == 1 or adjacency_matrix[k, l] == 100:

                        # Check if there is a cut between vertices i and j and vertices k and l
                        if (i <= k <= j < l) or (k < i <= l <= j):

                            # Record the cuts for each point
                            if index not in cut_huns:
                                cut_huns[index] = 1
                            else:
                                cut_huns[index] += 1

    # Print the cuts for each point
    for point, num_cuts in cut_huns.items():
        print(f"Point {point}: {num_cuts} cuts")
    cuthuns_array = np.array(list(cut_huns.values()))
    avg_error = sum(cuthuns_array) / len(cuthuns_array)
    # ............................plot

    cuts_ones = cutones_array
    cuts_hundreds = cuthuns_array

    # plotting the cuts for point ones
    # data = [cuts_ones, cuts_hundreds]
    data = [cuts_ones]

    # plotting the box plot horizontally
    plt.boxplot(data, vert=False)

    # setting the y-axis labels and legend
    plt.yticks([1, 2], ['correct', 'errors'])
    plt.ylabel('Points')
    plt.legend(['correct', 'errors'], loc='lower right')

    # setting the axis labels and title
    plt.xlabel('Cuts Number')
    plt.title('Box Plot of Cut Numbers')

    # displaying the plot
    plt.show()


    # Heatmap for deepbreaks prediction

    # Generate a random matrix with 1s and 0s
    matrix = pre

    # Define the list of points to highlight
    import itertools
    from matplotlib.colors import ListedColormap

    lst = [56, 203, 204, 205, 206, 207, 208, 425, 429, 430, 431, 467, 468, 469, 470, 471, 472]
    pairs = list(itertools.permutations(lst, 2))

    # Create a mask for the points to highlight
    mask = np.zeros_like(matrix)
    for point in pairs:
        mask[point[0], point[1]] = 1

    colors = ['#ffffff', '#f7f7f7', '#313695']
    cmap = ListedColormap(colors)

    # Create the heatmap plot
    fig, ax = plt.subplots()
    heatmap = ax.imshow(matrix, cmap='gray')
    highlight = ax.imshow(mask, cmap=cmap, alpha=0.5)
    ax.set_title('Heatmap with Highlighted Points')
    ax.axis('off')
    plt.show()


    # Get the coordinates of the 1s in the matrix
    coords = np.where(matrix == 1)
    coords = list(zip((coords[0]+1), (coords[1]+1)))

    # Print the coordinates
    print(f"The coordinates of the 1's in the matrix are: {coords}")

    included_pairs = []
    for pair in pairs:
        if pair in coords:
            included_pairs.append(pair)

    print(included_pairs)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass




#cut
z2 = z1
z2[z2 == 100] = 0
adjacency_matrix = z2
ones_indices = []
hundreds_indices = []

    # Iterate over the adjacency matrix and record the indices of the 1s and 100s
for i in range(adjacency_matrix.shape[0]):
        for j in range(adjacency_matrix.shape[1]):
            if adjacency_matrix[i, j] == 1:
                ones_indices.append((i, j))
            elif adjacency_matrix[i, j] == 100:
                hundreds_indices.append((i, j))

    # Initialize a dictionary to record the cuts for each point
cut_ones = {}
cut_huns = {}
    # Iterate over the ones and hundreds indices and check for cuts
for index in ones_indices:
        i, j = index

        # Check if there is an edge between vertices i and j
        if adjacency_matrix[i, j] == 1 or adjacency_matrix[i, j] == 100:

            # Iterate over all other pairs of vertices in the graph
            for k in range(adjacency_matrix.shape[0]):
                for l in range(k + 1, adjacency_matrix.shape[0]):

                    # Check if there is an edge between vertices k and l
                    if adjacency_matrix[k, l] == 1 or adjacency_matrix[k, l] == 100:

                        # Check if there is a cut between vertices i and j and vertices k and l
                        if (i <= k <= j < l) or (k < i <= l <= j):

                            # Record the cuts for each point
                            if index not in cut_ones:
                                cut_ones[index] = 1
                            else:
                                cut_ones[index] += 1

    # Print the cuts for each point
for point, num_cuts in cut_ones.items():
        print(f"Point {point}: {num_cuts} cuts")
cutones_array = np.array(list(cut_ones.values()))
avg_correct = sum(cutones_array) / len(cutones_array)
print(len(cutones_array))
print(avg_correct )

cuts_ones = cutones_array


    # plotting the cuts for point ones
    # data = [cuts_ones, cuts_hundreds]
data = [cuts_ones]

    # plotting the box plot horizontally
plt.boxplot(data, vert=False)

    # setting the y-axis labels and legend
plt.yticks([1, 2], ['correct', 'errors'])
plt.ylabel('Points')
plt.legend(['correct', 'errors'], loc='lower right')

    # setting the axis labels and title
plt.xlabel('Cuts Number')
plt.title('Box Plot of Cut Numbers')

    # displaying the plot
plt.show()