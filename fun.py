# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


# 密码子频率表
codon_frequencies = {
    'M': {'ATG': 1.0000},
    'V': {'GTG': 0.8571, 'GTA': 0.0714, 'GTC': 0.0714},
    'S': {'AGC': 0.0909, 'TCC': 0.9091},
    'K': {'AAG': 1.0000},
    'G': {'GGC': 0.9600, 'GGT': 0.0400},
    'E': {'GAG': 0.8571, 'GAA': 0.1429},
    'A': {'GCA': 0.1000, 'GCC': 0.7000, 'GCG': 0.2000},
    'I': {'ATC': 0.8571, 'ATT': 0.1429},
    'F': {'TTC': 1.0000},
    'R': {'CGG': 0.1818, 'CGC': 0.7273, 'AGG': 0.0909},
    'H': {'CAC': 1.0000},
    'N': {'AAC': 1.0000},
    'P': {'CCC': 0.7692, 'CCT': 0.2308},
    'Y': {'TAC': 0.9167, 'TAT': 0.0833},
    'T': {'ACC': 0.9375, 'ACA': 0.0625},
    'Q': {'CAG': 1.0000},
    'L': {'CTG': 0.7692, 'CTC': 0.0769, 'TTG': 0.1538},
    'W': {'TGG': 1.0000},
    'D': {'GAC': 1.0000}
}


# 遗传密码表（标准遗传密码）
genetic_code = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
}

# 反转遗传密码表，以便从氨基酸到密码子
reverse_genetic_code = {}
for codon, aa in genetic_code.items():
    if aa not in reverse_genetic_code:
        reverse_genetic_code[aa] = []
    reverse_genetic_code[aa].append(codon)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    dna_sequence = input("请输入DNA序列：")

    # 原始DNA序列
    # dna_sequence = "atggccccaaagaagaagcggaaggtcggtatccacggagtcccagcagccaagcggaactacatcctgggcctggacatcggcatcaccagcgtgggctacggcatcatcgactacgagacacgggacgtgatcgatgccggcgtgcggctgttcaaagaggccaacgtggaaaacaacgagggcaggcggagcaagagaggcgccagaaggctgaagcggcggaggcggcatagaatccagagagtgaagaagctgctgttcgactacaacctgctgaccgaccacagcgagctgagcggcatcaacccctacgaggccagagtgaagggcctgagccagaagctgagcgaggaagagttctctgccgccctgctgcacctggccaagagaagaggcgtgcacaacgtgaacgaggtggaagaggacaccggcaacgagctgtccaccaaagagcagatcagccggaacagcaaggccctggaagagaaatacgtggccgaactgcagctggaacggctgaagaaagacggcgaagtgcggggcagcatcaacagattcaagaccagcgactacgtgaaagaagccaaacagctgctgaaggtgcagaaggcctaccaccagctggaccagagcttcatcgacacctacatcgacctgctggaaacccggcggacctactatgagggacctggcgagggcagccccttcggctggaaggacatcaaagaatggtacgagatgctgatgggccactgcacctacttccccgaggaactgcggagcgtgaagtacgcctacaacgccgacctgtacaacgccctgaacgacctgaacaatctcgtgatcaccagggacgagaacgagaagctggaatattacgagaagttccagatcatcgagaacgtgttcaagcagaagaagaagcccaccctgaagcagatcgccaaagaaatcctcgtgaacgaagaggatattaagggctacagagtgaccagcaccggcaagcccgagttcaccaacctgaaggtgtaccacgacatcaaggacattaccgcccggaaagagattattgagaacgccgagctgctggatcagattgccaagatcctgaccatctaccagagcagcgaggacatccaggaagaactgaccaatctgaactccgagctgacccaggaagagatcgagcagatctctaatctgaagggctataccggcacccacaacctgagcctgaaggccatcaacctgatcctggacgagctgtggcacaccaacgacaaccagatcgctatcttcaaccggctgaagctggtgcccaagaaggtggacctgtcccagcagaaagagatccccaccaccctggtggacgacttcatcctgagccccgtcgtgaagagaagcttcatccagagcatcaaagtgatcaacgccatcatcaagaagtacggcctgcccaacgacatcattatcgagctggcccgcgagaagaactccaaggacgcccagaaaatgatcaacgagatgcagaagcggaaccggcagaccaacgagcggatcgaggaaatcatccggaccaccggcaaagagaacgccaagtacctgatcgagaagatcaagctgcacgacatgcaggaaggcaagtgcctgtacagcctggaagccatccctctggaagatctgctgaacaaccccttcaactatgaggtggaccacatcatccccagaagcgtgtccttcgacaacagcttcaacaacaaggtgctcgtgaagcaggaagaaaacagcaagaagggcaaccggaccccattccagtacctgagcagcagcgacagcaagatcagctacgaaaccttcaagaagcacatcctgaatctggccaagggcaagggcagaatcagcaagaccaagaaagagtatctgctggaagaacgggacatcaacaggttctccgtgcagaaagacttcatcaaccggaacctggtggataccagatacgccaccagaggcctgatgaacctgctgcggagctacttcagagtgaacaacctggacgtgaaagtgaagtccatcaatggcggcttcaccagctttctgcggcggaagtggaagtttaagaaagagcggaacaaggggtacaagcaccacgccgaggacgccctgatcattgccaacgccgatttcatcttcaaagagtggaagaaactggacaaggccaaaaaagtgatggaaaaccagatgttcgaggaaaagcaggccgagagcatgcccgagatcgaaaccgagcaggagtacaaagagatcttcatcaccccccaccagatcaagcacattaaggacttcaaggactacaagtacagccaccgggtggacaagaagcctaatagagagctgattaacgacaccctgtactccacccggaaggacgacaagggcaacaccctgatcgtgaacaatctgaacggcctgtacgacaaggacaatgacaagctgaaaaagctgatcaacaagagccccgaaaagctgctgatgtaccaccacgacccccagacctaccagaaactgaagctgattatggaacagtacggcgacgagaagaatcccctgtacaagtactacgaggaaaccgggaactacctgaccaagtactccaaaaaggacaacggccccgtgatcaagaagattaagtattacggcaacaaactgaacgcccatctggacatcaccgacgactaccccaacagcagaaacaaggtcgtgaagctgtccctgaagccctacagattcgacgtgtacctggacaatggcgtgtacaagttcgtgaccgtgaagaatctggatgtgatcaaaaaagaaaactactacgaagtgaatagcaagtgctatgaggaagctaagaagctgaagaagatcagcaaccaggccgagtttatcgcctccttctacaacaacgatctgatcaagatcaacggcgagctgtatagagtgatcggcgtgaacaacgacctgctgaaccggatcgaagtgaacatgatcgacatcacctaccgcgagtacctggaaaacatgaacgacaagaggccccccaggatcattaagacaatcgcctccaagacccagagcattaagaagtacagcacagacattctgggcaacctgtatgaagtgaaatctaagaagcaccctcagatcatcaaaaagggcaaaaggccggcggccacgaaaaaggccggccaggcaaaaaagaaaaagggatcctacccatacgatgttccagattacgcttacccatacgatgttccagattacgcttacccatacgatgttccagattacgcttaa"
    # 确保序列长度是3的倍数
    if len(dna_sequence) % 3 != 0:
        raise ValueError("The length of the DNA sequence must be a multiple of 3.")

    # 分割成密码子
    original_codons = [dna_sequence[i:i+3] for i in range(0, len(dna_sequence), 3)]

    # 将原始DNA序列翻译成氨基酸序列
    amino_acids = [genetic_code[codon.upper()] for codon in original_codons]

    # 优化后的密码子列表
    optimized_codons = []

    # 选择最优的密码子
    for codon, aa in zip(original_codons, amino_acids):
        if aa in codon_frequencies:
            # 选择频率最高的密码子
            best_codon = max(codon_frequencies[aa], key=codon_frequencies[aa].get)
            optimized_codons.append(best_codon)
        else:
            # 如果没有提供该氨基酸的频率信息，保留原始密码子
            optimized_codons.append(codon.upper())

    # 重新构建优化后的DNA序列
    optimized_dna_sequence = ''.join(optimized_codons)

    # 输出结果
    print("Original DNA Sequence:")
    print(dna_sequence)
    print("\nOptimized DNA Sequence:")
    print(optimized_dna_sequence)