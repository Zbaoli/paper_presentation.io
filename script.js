document.addEventListener('DOMContentLoaded', () => {
    const codeInput = document.getElementById('code-input');
    const calculateBtn = document.getElementById('calculate-btn');
    const outputDisplay = document.getElementById('output-display');
    
    calculateBtn.addEventListener('click', () => {
        try {
            // 获取输入的代码
            const code = codeInput.value.trim();
            
            if (!code) {
                outputDisplay.textContent = '请输入代码';
                return;
            }
            
            // 设置一个安全的计算环境
            const sandbox = {
                result: null,
                console: {
                    log: function(...args) {
                        sandbox.result = (sandbox.result || '') + args.join(' ') + '\n';
                    }
                }
            };
            
            // 包装用户代码，捕获计算结果
            const wrappedCode = `
                try {
                    ${code}
                } catch (error) {
                    console.log("计算错误:", error.message);
                }
            `;
            
            // 使用Function构造器创建一个新函数
            const computeFunction = new Function(...Object.keys(sandbox), wrappedCode);
            
            // 执行计算
            computeFunction(...Object.values(sandbox));
            
            // 显示结果
            if (sandbox.result !== null) {
                outputDisplay.textContent = sandbox.result;
            } else {
                // 如果代码没有使用console.log，尝试直接评估表达式
                try {
                    const directResult = eval(code);
                    outputDisplay.textContent = directResult !== undefined ? directResult : '计算完成，但没有返回值';
                } catch (error) {
                    outputDisplay.textContent = '计算错误: ' + error.message;
                }
            }
        } catch (error) {
            outputDisplay.textContent = '系统错误: ' + error.message;
        }
    });
    
    // 添加示例代码
    codeInput.value = `// 示例：计算斐波那契数列
function fibonacci(n) {
    if (n <= 1) return n;
    return fibonacci(n-1) + fibonacci(n-2);
}

// 计算前10个斐波那契数
for (let i = 0; i < 10; i++) {
    console.log(\`fibonacci(\${i}) = \${fibonacci(i)}\`);
}`;
});