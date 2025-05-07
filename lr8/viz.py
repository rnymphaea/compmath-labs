import sys
import matplotlib.pyplot as plt
import re

def parse_data(input_lines):
    data = []
    # Шаблон для парсинга чисел с научной нотацией
    pattern = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
    
    for line in input_lines:
        # Пропускаем разделительные строки
        if '----------' in line or not line.strip():
            continue
        
        # Ищем все числа в строке
        matches = re.findall(pattern, line)
        if len(matches) >= 5:
            try:
                x = float(matches[0])
                eq = float(matches[1])
                ch = float(matches[2])
                orig = float(matches[3])
                data.append((x, eq, ch, orig))
            except ValueError:
                continue
    return data

def plot_data(data):
    x = [d[0] for d in data]
    eq = [d[1] for d in data]
    ch = [d[2] for d in data]
    orig = [d[3] for d in data]

    plt.figure(figsize=(12, 8))
    
    # График функций
    plt.subplot(2, 1, 1)
    plt.plot(x, orig, 'k-', label='Исходная функция', linewidth=2)
    plt.plot(x, eq, 'r--', label='Равноотстоящие узлы')
    plt.plot(x, ch, 'b-.', label='Узлы Чебышева')
    plt.title('Сравнение интерполяции')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.grid(True)
    plt.legend()

    # График отклонений
    plt.subplot(2, 1, 2)
    plt.plot(x, [eq[i] - orig[i] for i in range(len(x))], 'r--', label='Отклонение (равноот.)')
    plt.plot(x, [ch[i] - orig[i] for i in range(len(x))], 'b-.', label='Отклонение (Чебышев)')
    plt.title('Отклонения интерполяции')
    plt.xlabel('x')
    plt.ylabel('Δf(x)')
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    print("Введите данные (Ctrl+D для завершения ввода):")
    input_lines = sys.stdin.read().splitlines()
    data = parse_data(input_lines)
    
    if not data:
        print("Не удалось распознать данные")
        sys.exit(1)
    
    plot_data(data)
