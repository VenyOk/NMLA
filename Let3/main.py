import matplotlib.pyplot as plt
import numpy as np

# Данные для варианта 31
frequencies = [63, 125, 250, 500, 1000, 2000, 4000, 8000]
L_calculated = [67.19, 71.06, 69.37, 71.32, 69.78, 67.28, 76.48, 78.35]
L_normative = [83, 74, 68, 63, 60, 57, 55, 54]

# Расчет требуемого снижения шума
required_reduction = [max(0, calc - norm) for calc, norm in zip(L_calculated, L_normative)]

print("Требуемое снижение уровней шума по частотам:")
for freq, reduction in zip(frequencies, required_reduction):
    print(f"{freq} Гц: {reduction:.1f} дБ")

# Создаем фигуру с двумя графиками
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))

# График 1: Сравнение спектров
ax1.semilogx(frequencies, L_calculated, 'bo-', linewidth=3, markersize=10, label='Расчетный спектр (Lрасч)')
ax1.semilogx(frequencies, L_normative, 'ro-', linewidth=3, markersize=10, label='Предельный спектр (Lдоп)')

# Настройка первого графика
ax1.set_xlabel('Частота, Гц', fontsize=12)
ax1.set_ylabel('Уровень звукового давления, дБ', fontsize=12)
ax1.set_title('Сравнение нормативных и расчетных уровней звукового давления\nВариант 31', fontsize=14, fontweight='bold')
ax1.grid(True, which='both', linestyle='--', alpha=0.7)
ax1.legend(fontsize=12)
ax1.set_xlim(50, 10000)
ax1.set_ylim(50, 90)

# Добавляем подписи значений
for i, (f, l_calc, l_norm) in enumerate(zip(frequencies, L_calculated, L_normative)):
    ax1.annotate(f'{l_calc}', (f, l_calc), textcoords="offset points", xytext=(0,10), 
                ha='center', fontsize=10, fontweight='bold', color='blue')
    ax1.annotate(f'{l_norm}', (f, l_norm), textcoords="offset points", xytext=(0,-15), 
                ha='center', fontsize=10, fontweight='bold', color='red')

# Выделяем зону превышения
for i in range(len(frequencies)):
    if L_calculated[i] > L_normative[i]:
        ax1.fill_between([frequencies[i]*0.9, frequencies[i]*1.1], 
                        [L_calculated[i], L_calculated[i]], 
                        [L_normative[i], L_normative[i]], 
                        alpha=0.3, color='red')

# График 2: Требуемое снижение шума
bars = ax2.bar(range(len(frequencies)), required_reduction, color='orange', alpha=0.7)
ax2.set_xlabel('Частота, Гц', fontsize=12)
ax2.set_ylabel('Требуемое снижение ΔL, дБ', fontsize=12)
ax2.set_title('Требуемое снижение уровней шума', fontsize=14, fontweight='bold')
ax2.set_xticks(range(len(frequencies)))
ax2.set_xticklabels(frequencies)
ax2.grid(True, axis='y', linestyle='--', alpha=0.7)

# Добавляем значения на столбцы
for bar, reduction in zip(bars, required_reduction):
    if reduction > 0:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{reduction:.1f} дБ', ha='center', va='bottom', fontweight='bold')

plt.tight_layout()
plt.savefig('noise_reduction_analysis_variant31.png', dpi=300, bbox_inches='tight')
plt.show()

# Дополнительный детализированный график
plt.figure(figsize=(12, 8))

# Основной график сравнения
plt.semilogx(frequencies, L_calculated, 'bo-', linewidth=3, markersize=10, label='Lрасч (расчетный)')
plt.semilogx(frequencies, L_normative, 'ro-', linewidth=3, markersize=10, label='Lдоп (нормативный)')

# Заполнение области превышения
for i in range(len(frequencies)-1):
    if L_calculated[i] > L_normative[i] or L_calculated[i+1] > L_normative[i+1]:
        x_fill = [frequencies[i], frequencies[i+1]]
        y1_fill = [L_calculated[i], L_calculated[i+1]]
        y2_fill = [L_normative[i], L_normative[i+1]]
        plt.fill_between(x_fill, y1_fill, y2_fill, alpha=0.3, color='red', label='Превышение' if i == 3 else "")

plt.xlabel('Частота, Гц', fontsize=12)
plt.ylabel('Уровень звукового давления, дБ', fontsize=12)
plt.title('Сравнение нормативных и расчетных уровней звукового давления\nВариант 31', fontsize=14, fontweight='bold')
plt.grid(True, which='both', linestyle='--', alpha=0.7)
plt.legend(fontsize=12)
plt.xlim(50, 10000)
plt.ylim(50, 90)

# Добавляем аннотации с превышениями
for i, (f, calc, norm) in enumerate(zip(frequencies, L_calculated, L_normative)):
    if calc > norm:
        plt.annotate(f'ΔL = {calc-norm:.1f} дБ', (f, (calc+norm)/2), 
                    textcoords="offset points", xytext=(0,0), 
                    ha='center', fontsize=10, fontweight='bold',
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7))

plt.tight_layout()
plt.savefig('detailed_noise_comparison_variant31.png', dpi=300, bbox_inches='tight')
plt.show()

# Создаем таблицу требуемого снижения
print("\n" + "="*60)
print("ТАБЛИЦА ТРЕБУЕМОГО СНИЖЕНИЯ УРОВНЕЙ ШУМА")
print("="*60)
print("Частота, Гц | Расчетный уровень | Нормативный уровень | Требуемое снижение")
print("-"*60)
for freq, calc, norm in zip(frequencies, L_calculated, L_normative):
    reduction = max(0, calc - norm)
    print(f"{freq:^11} | {calc:^16.1f} | {norm:^18.1f} | {reduction:^16.1f}")
print("="*60)