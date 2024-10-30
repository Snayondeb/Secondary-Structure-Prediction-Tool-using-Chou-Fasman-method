# On the basis of Propensity Value(Chou-Fasman Method) prediction of Protein Secondary Structure
import numpy as np
import requests
import plotext as plt
from tabulate import tabulate
import sys
from colorama import Fore, Style

#Creating APIs
WEBSITE_API = "https://rest.uniprot.org/"

PROTEINS_API = "https://www.ebi.ac.uk/proteins/api"

print("On the basis of Propensity Value(Chou-Fasman Method) prediction of Protein Secondary Structure")

# Amino Acid Abbreviations
amino_acids = {
    'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine',
    'D': 'Aspartic Acid', 'C': 'Cysteine', 'E': 'Glutamic Acid',
    'Q': 'Glutamine', 'G': 'Glycine', 'H': 'Histidine',
    'I': 'Isoleucine', 'L': 'Leucine', 'K': 'Lysine',
    'M': 'Methionine', 'F': 'Phenylalanine', 'P': 'Proline',
    'S': 'Serine', 'T': 'Threonine', 'W': 'Tryptophan',
    'Y': 'Tyrosine', 'V': 'Valine'
}

# Chou-Fasman Parameters
parameters = {
    'Alpha Helix': {'A': 1.42, 'R': 0.98, 'N': 0.67, 'D': 0.73, 'C': 0.70, 'E': 1.51, 'Q': 1.11, 'G': 0.57, 'H': 1.00, 'I': 1.08, 'L': 1.21, 'K': 1.16, 'M': 1.45, 'F': 1.13, 'P': 0.57, 'S': 0.76, 'T': 0.82, 'W': 1.08, 'Y': 0.69, 'V': 1.06},
    'Beta Sheet': {'A': 0.83, 'R': 0.93, 'N': 0.89, 'D': 1.01, 'C': 1.19, 'E': 0.74, 'Q': 0.80, 'G': 0.81, 'H': 0.81, 'I': 1.60, 'L': 1.22, 'K': 0.81, 'M': 1.05, 'F': 1.28, 'P': 0.62, 'S': 1.13, 'T': 1.20, 'W': 1.14, 'Y': 1.47, 'V': 1.65},
    'Coil': {'A': 0.75, 'R': 1.09, 'N': 1.44, 'D': 1.26, 'C': 0.11, 'E': 0.75, 'Q': 1.09, 'G': 1.62, 'H': 1.19, 'I': 0.32, 'L': 0.57, 'K': 1.03, 'M': 0.50, 'F': 0.59, 'P': 1.81, 'S': 0.11, 'T': 0.98, 'W': 0.78, 'Y': 0.84, 'V': 0.29},
    'Turn' : {'A': 0.66, 'R': 0.95, 'N': 1.56, 'D': 1.46, 'C': 1.19, 'E': 0.74, 'Q': 0.98, 'G': 1.56, 'H': 0.95, 'I': 0.47, 'L': 0.59, 'K': 1.01, 'M': 0.60, 'F': 0.60, 'P': 1.52, 'S': 1.43, 'T': 0.96, 'W': 0.96, 'Y': 1.14, 'V': 0.50}
}

#function to Predict the Secondary Structure on the basis of Chou Fasman method
def calculate_structure(probability_matrix, sequence):
    structure_prediction = []
    for amino_acid in sequence:
        probabilities = [parameters['Alpha Helix'][amino_acid], parameters['Beta Sheet'][amino_acid],parameters['Coil'][amino_acid], parameters['Turn'][amino_acid]]
        max_index = probabilities.index(max(probabilities)) #probabilities = [0.98, 0.93]
        if max_index == 0:
            structure_prediction.append('H')  # Alpha Helix
        elif max_index == 1:
            structure_prediction.append('E')  # Beta Sheet
        elif max_index == 3:
            structure_prediction.append('T')  # Turn
        else:
            structure_prediction.append('C')  # Coil
    return ''.join(structure_prediction)

#main function that control all the task except prediction
def main():
    #sequence = input("Enter protein sequence: ")
    #'''
    print("1. If you know/have the FASTA sequence, go for the option 1\n2. If you don't know/have the FASTA sequence, go for the option 2")
    ch = int(input("Enter the choice for opting the way to predict the secondary structure of protein: "))
    if ch == 1:
        print("Paste the FASTA sequence. Press ENTER, then Ctrl+D to save it.")
        sequence = ''
        while True:
            try:
                line = input()
            except EOFError:
                break
            sequence += line

    elif ch == 2:
        def get_url(url, **kwargs):
            response = requests.get(url, **kwargs);
            if not response.ok:
                print(response.text)
                response.raise_for_status()
                sys.exit()
            return response

        protein_name = input("Enter the name of the protein: ")
        print("The types of taxonomy")
        print(
            "1. BOVIN for Bovine\n2. CHICK for Chicken\n3. ECOLI for Escherichia coli\n4. HORSE for Horse\n5. HUMAN for Homo sapiens\n6. MAIZE for Maize (Zea mays)\n7. MOUSE for Mouse\n8. PEA for Garden pea (Pisum sativum)\n9. PIG for Pig\n10. RABIT for Rabbit\n11. RAT for Rat\n12. SHEEP for Sheep\n13. SOYBN for Soybean (Glycine max)\n14. TOBAC for Common tobacco (Nicotina tabacum)\n15. WHEAT for Wheat (Triticum aestivum)\n16. YEAST for Baker's yeast (Saccharomyces cerevisiae)\n17. * for ALL(It will have all types of structure)")
        taxonomy_type = input("Enter the taxonomy of whose secondary protein structure is required: ").upper()
        r = get_url(
            f"{WEBSITE_API}/uniprotkb/search?query={protein_name} AND (taxonomy_name={taxonomy_type})&format=fasta")

        # creating a temporary text file and deleting it after the runtime.
        fhw = open("temp.txt", "w+")
        file = str(r.text)
        fhw.writelines(file)
        fhw.close()
        fhr = open("temp.txt", "r")
        st = ''
        lst_k = []
        lst_v = []
        # lst_check = [">sp|P01308|INS_HUMAN Insulin OS=Homo sapiens OX=9606 GN=INS PE=1 SV=1\n","MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED\n","LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN\n", ">sp|P06213|INSR_HUMAN Insulin receptor OS=Homo sapiens OX=9606 GN=INSR PE=1 SV=4\n", "MATGGRRGAAAAPLLVAVAALLLGAAGHLYPGEVCPGMDIRNNLTRLHELENCSVIEGHL\n", "QILLMFKTRPEDFRDLSFPKLIMITDYLLLFRVYGLESLKDLFPNLTVIRGSRLFFNYAL\n", "FEDMENVPLDRSSHCQREEAGGRDGGSSLGFKRSYEEHIPYTHMNGGKKNGRILTLPRSN\n","PS\n",">sp|P14735|IDE_HUMAN Insulin-degrading enzyme OS=Homo sapiens OX=9606 GN=IDE PE=1 SV=4\n","MRYRLAWLLHPALPSTFRSVLGARLPPPERLCGFQKKTYSKMNNPAIKRIGNHITKSPED\n"]
        for k in fhr.readlines():
            if k.startswith(">"):
                lst_v.append(st)
                st = ''
                lst_k.append(k.rstrip("\n"))
            elif k.isupper():
                st += k.rstrip("\n")
            else:
                pass
        lst_v.append(st)
        lst_v.pop(0)
        f_lst_k = f_lst_v = []
        flag = 0
        n_dic = {}
        dic = dict(zip(lst_k, lst_v))
        keys = list(dic.keys())
        values = []  # list(dic.values())
        for k in range(0, len(lst_k)):
            key = lst_k[k]
            values.append(int((((key.split(" "))[-2]).split("="))[-1]))
        arr = np.array(values)
        sorted_value_index = (-arr).argsort()
        rev_sorted_dict = {keys[i]: values[i] for i in sorted_value_index}
        #print(rev_sorted_dict)
        n_keys = list(rev_sorted_dict.keys())
        # print(n_keys)
        for j in n_keys:
            n_dic[j] = dic[j]
        # print(n_dic)
        counter_1 = 0
        for head in n_dic:
            print(Fore.YELLOW + "(",counter_1 + 1, ") ", head)
            print(Style.RESET_ALL, end='')
            print(n_dic[head])
            counter_1 += 1
        print(Fore.LIGHTRED_EX+"Note:- More the PE values, more is the Existence of the Protein\n\t   The list has been ordered as Highest PE value at the top, and decreases as we go down.")
        print(Style.RESET_ALL, end='')
        ch_2 = int(input("Enter the serial number of the protein: "))
        sequence = n_dic[n_keys[ch_2-1]]
        fhr.close()


    sequence = sequence.upper()
    print("No. of Amino acids : ",len(sequence))
    prediction = calculate_structure(parameters, sequence)
    print(Fore.CYAN + "Predicted Secondary Structure: ", prediction)

    t = len(prediction)
    count_H, count_E, count_C, count_T = 0, 0, 0, 0
    for k in prediction:
        if k == "H":
            count_H += 1
        elif k == "E":
            count_E += 1
        elif k == "C":
            count_C += 1
        elif k == "T":
            count_T += 1
        else:
            pass
    percent_H = (count_H / t) * 100
    percent_E = (count_E / t) * 100
    percent_C = (count_C / t) * 100
    percent_T = (count_T / t) * 100
    #print(percent_H, "\n", percent_E, "\n", percent_C, "\n", percent_T)
    data = [["alpa-helix(H)", percent_H], ["beta-sheet(E)", percent_E], ["coil(C)", percent_C], ["turn(T)", percent_T]]
    print(Fore.YELLOW + tabulate(data, headers=["Types of Bond", "Percentage"]))
    x_axis = ["alpa-helix(H)", "beta-sheet(E)", "coil(C)", "turn(T)"]
    y_axis = [percent_H, percent_E, percent_C, percent_T]

    plt.simple_bar(x_axis, y_axis, width=100, title='The Graphical representation of % of helix, sheet, turn, coil')
    plt.show()
    print(Fore.LIGHTRED_EX+"\t1. Red - Alpha-Helix(H)")
    print(Fore.LIGHTGREEN_EX + "\t2. Green - Beta-Sheet(E)")
    print(Fore.LIGHTBLUE_EX + "\t3. Blue - Coil(C)")
    print(Fore.LIGHTYELLOW_EX + "\t4. Yellow - Terminal(T)")
    print(Style.BRIGHT+Fore.LIGHTWHITE_EX+"\tTherefore, the Actual Secondary Structure is:-")
    counter = 0
    for k in prediction:
        if k == "H":
            print(Fore.LIGHTRED_EX + f'{sequence[counter]}', end='')
        elif k == "E":
            print(Fore.LIGHTGREEN_EX + f'{sequence[counter]}', end='')
        elif k == "C":
            print(Fore.LIGHTBLUE_EX + f'{sequence[counter]}', end='')
        elif k == "T":
            print(Fore.LIGHTYELLOW_EX + f'{sequence[counter]}', end='')
        else:
            pass
        counter += 1
        if counter%100 == 0:
            print("\n",end='')


if __name__ == "__main__":
    main()
