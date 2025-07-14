// 语言切换功能
document.addEventListener('DOMContentLoaded', function() {
    const langSwitchBtn = document.getElementById('lang-switch');
    let currentLang = 'en';

    // 设置输入框的占位符文本
    const codeInput = document.getElementById('code-input');
    updatePlaceholder();

    // 语言切换按钮点击事件
    langSwitchBtn.addEventListener('click', function() {
        currentLang = currentLang === 'en' ? 'zh' : 'en';
        switchLanguage();
    });

    function switchLanguage() {
        // 切换所有带有data-lang属性的元素
        document.querySelectorAll('[data-lang]').forEach(element => {
            if (element.getAttribute('data-lang') === currentLang) {
                element.style.display = '';
            } else {
                element.style.display = 'none';
            }
        });

        // 更新输入框占位符
        updatePlaceholder();
    }

    function updatePlaceholder() {
        codeInput.placeholder = currentLang === 'en' 
            ? 'Please enter DNA sequence here...' 
            : '请输入DNA序列...';
    }

    // 优化按钮点击事件
    document.getElementById('calculate-btn').addEventListener('click', function() {
        const input = codeInput.value.trim();
        if (!input) {
            alert(currentLang === 'en' ? 'Please enter a DNA sequence.' : '请输入DNA序列。');
            return;
        }

        // 这里添加你的优化逻辑
        const output = `Optimized sequence for ${input}`;
        document.getElementById('output-display').textContent = output;
    });
});

// 密码子频率表
const codonFrequencies = {
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
};

// 遗传密码表（标准遗传密码）
const geneticCode = {
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
};

// 反转遗传密码表，以便从氨基酸到密码子
const reverseGeneticCode = {};
for (const [codon, aa] of Object.entries(geneticCode)) {
    if (!reverseGeneticCode[aa]) {
        reverseGeneticCode[aa] = [];
    }
    reverseGeneticCode[aa].push(codon);
}

// DNA序列优化函数
function optimizeDNA(dnaSequence) {
    // 确保序列长度是3的倍数
    if (dnaSequence.length % 3 !== 0) {
        throw new Error("DNA sequence length must be a multiple of 3");
    }

    // 分割成密码子
    const originalCodons = [];
    for (let i = 0; i < dnaSequence.length; i += 3) {
        originalCodons.push(dnaSequence.slice(i, i + 3));
    }

    // 将原始DNA序列翻译成氨基酸序列
    const aminoAcids = originalCodons.map(codon => geneticCode[codon.toUpperCase()]);

    // 优化后的密码子列表
    const optimizedCodons = [];

    // 选择最优的密码子
    for (let i = 0; i < originalCodons.length; i++) {
        const codon = originalCodons[i];
        const aa = aminoAcids[i];
        
        if (codonFrequencies[aa]) {
            // 选择频率最高的密码子
            const bestCodon = Object.entries(codonFrequencies[aa])
                .reduce((a, b) => a[1] > b[1] ? a : b)[0];
            optimizedCodons.push(bestCodon);
        } else {
            // 如果没有提供该氨基酸的频率信息，保留原始密码子
            optimizedCodons.push(codon.toUpperCase());
        }
    }

    // 重新构建优化后的DNA序列
    return optimizedCodons.join('');
}

// 页面加载完成后执行
document.addEventListener('DOMContentLoaded', () => {
    const codeInput = document.getElementById('code-input');
    const calculateBtn = document.querySelector('#calculate-btn[data-lang="en"]');
    const outputDisplay = document.getElementById('output-display');
    
    // 添加示例DNA序列
    codeInput.value = ``;
    
    calculateBtn.addEventListener('click', () => {
        try {
            const dnaSequence = codeInput.value.trim();
            
            if (!dnaSequence) {
                outputDisplay.textContent = 'Please enter a DNA sequence';
                return;
            }
            
            // 优化DNA序列
            const optimizedSequence = optimizeDNA(dnaSequence);
            
            // 显示结果
            // outputDisplay.textContent = `Original DNA Sequence:\n${dnaSequence}\n\nOptimized DNA Sequence:\n${optimizedSequence}`;
            outputDisplay.textContent = `${optimizedSequence}`;
        } catch (error) {
            outputDisplay.textContent = 'Error: ' + error.message;
        }
    });
});