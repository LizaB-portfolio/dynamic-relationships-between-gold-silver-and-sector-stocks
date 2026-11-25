import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels
from itertools import combinations
from scipy import stats
from pathlib import Path
import statsmodels.api as sm

from statsmodels.stats.multitest import multipletests

# Очистка и подготовка данных
FILES = {
    'gold.xlsx': 'Gold',
    'silver.xlsx': 'Silver',
    'Consumer_discretionary.xlsx': 'Cons_discr',
    'Consumer_staples.xlsx': 'Cons_staples',
    'Energy.xlsx': 'Energy',
    'Financials.xlsx': 'Financials',
    'Information_technology.xlsx': 'Inf_tech',
    'Health.xlsx': 'Health',
    'Industrials.xlsx': 'Industrials',
    'Materials.xlsx': 'Materials',
    'Real_estate.xlsx': 'Real_estate',
    'Telecom.xlsx': 'Telecom',
    'Utilities.xlsx': 'Utilities',
}

GROUPS = {
    'c': ['Gold', 'Silver', 'Cons_discr', 'Financials', 'Industrials', 'Materials'],  # циклические
    'd': ['Gold', 'Silver', 'Cons_staples', 'Health', 'Utilities'],                   # защитные
    'm': ['Gold', 'Silver', 'Inf_tech', 'Real_estate', 'Telecom', 'Energy'],          # смешанные
}

# Путь к папке с исходными Excel-файлами
BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / 'data_raw'

def load_series(file: str, col: str) -> pd.DataFrame:
    file_path = DATA_DIR / file
    df = pd.read_excel(file_path)
    df.columns = ['Date', col]
    df['Date'] = pd.to_datetime(df['Date'])
    df[col] = pd.to_numeric(
        df[col].astype(str).str.replace(',', '.'),
        errors='coerce'
    )
    return df.set_index('Date')

series = {col: load_series(file, col) for file, col in FILES.items()}
panels = {}
for suffix, cols in GROUPS.items():
    panel = pd.concat([series[c] for c in cols], axis=1).sort_index().dropna()
    assert panel.index.is_monotonic_increasing, 'Даты не в хронологическом порядке!'
    panels[suffix] = panel

final_levels = {s: np.log(p) for s, p in panels.items()}
final_returns = {s: (final_levels[s].diff() * 100).dropna() for s in panels}
model3 = pd.concat(series.values(), axis=1).sort_index()
model3_c, model3_d, model3_m = panels['c'], panels['d'], panels['m']
final_levels_c, final_levels_d, final_levels_m = final_levels['c'], final_levels['d'], final_levels['m']
final_returns_c, final_returns_d, final_returns_m = final_returns['c'], final_returns['d'], final_returns['m']

for s in ('c', 'd', 'm'):
    print(f'Группа {s}: {panels[s].shape[0]} наблюдений, '
          f'{panels[s].index[0].date()} - {panels[s].index[-1].date()}')

# EDA
RETURNS = {'cyclical': final_returns_c,
           'defensive': final_returns_d,
           'mixed': final_returns_m}

OR_PRICES = {'cyclical': model3_c, 'defensive': model3_d, 'mixed': model3_m}

# Путь к папке для графиков
FIGURES_DIR = BASE_DIR / 'outputs' / 'figures'

# Автоматическое создание папки, если она отсутствует
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Визуализация динамики цен
# График нормированных цен с исключением информационных технологий (для таблицы смешанных секторов, для циклических и защитных исключенных столбцов нет)
for group_name, df in OR_PRICES.items():

    # Исключаем Inf_tech только для смешанной группы
    if group_name == 'mixed':
        data = df.drop(columns=['Inf_tech'])
    else:
        data = df.copy()

    # Нормирование цен
    model_indexed = data.div(data.iloc[0]) * 100

    plt.figure(figsize=(14, 7))

    for col in model_indexed.columns:
        plt.plot(model_indexed.index, model_indexed[col], label=col)

    plt.axhline(y=100, color='black', linestyle='--', linewidth=0.8)

    plt.legend(loc='best')
    plt.title(f'Динамика нормированных цен активов: {group_name} (начальное значение = 100)')
    plt.xlabel('Дата')
    plt.ylabel('Индекс')
    plt.grid(alpha=0.3)

    # Сохранение графика
    plt.savefig(FIGURES_DIR / f'normalized_prices_{group_name}.png',
                dpi=300, bbox_inches='tight')

    plt.show()
    plt.close()

# Обычный график, показывающий динамику цен (для отрасли информационных технологий)
plt.figure(figsize=(12,6))
plt.plot(model3_m['Inf_tech'])
 
plt.title('Информационные технологии, долл. США')
plt.xlabel('Дата')
plt.ylabel('Inf_tech')
plt.grid(True)

plt.savefig(FIGURES_DIR / 'Inf_tech_prices.png',
            dpi=300, bbox_inches='tight')

plt.show()
plt.close()

# Графики доходности (разница логарифмов)
for group_name, df in RETURNS.items():
    for col in df.columns:

        plt.figure(figsize=(12, 6))
        plt.plot(df[col])

        plt.title(f'{col}, % ({group_name})')
        plt.xlabel('Дата')
        plt.ylabel(col)
        plt.grid(True)

        filename = f'{col}_{group_name}_returns.png'
        plt.savefig(FIGURES_DIR / filename,
                    dpi=300, bbox_inches='tight')

        plt.show()
        plt.close()

# Описательная статистика
def descriptive_stats(returns: pd.DataFrame) -> pd.DataFrame:
    out = returns.describe().T
    out['skew'] = returns.skew()
    out['kurtosis'] = returns.kurtosis()
    out['jarque_bera_p'] = returns.apply(lambda s: stats.jarque_bera(s)[1])
    return out.T


pd.set_option('display.width', None)
pd.set_option('display.max_columns', None)

for name, df in RETURNS.items():

    stats_df = descriptive_stats(df).round(4)

    print(f'\n===== Описательная статистика, группа {name} =====')
    print(stats_df)

    fig, ax = plt.subplots(figsize=(14, 7))
    ax.axis('off')

    table = ax.table(
        cellText=stats_df.values,
        rowLabels=stats_df.index,
        colLabels=stats_df.columns,
        cellLoc='center',
        loc='center'
    )

    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.2, 1.5)

    plt.title(f'Описательная статистика: {name}', pad=20)

    filename = f'descriptive_stats_{name}.png'
    plt.savefig(
        FIGURES_DIR / filename,
        dpi=300,
        bbox_inches='tight'
    )

    plt.show()
    plt.close()

# Корреляционные матрицы
def corr_with_pvalues(
    df: pd.DataFrame,
    hac_maxlags: int = 5
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Рассчитывает:
    1. Матрицу коэффициентов корреляции Пирсона.
    2. Матрицу необработанных HAC p-value.
    3. Матрицу p-value с поправкой Холма.

    Для проверки H0: отсутствие линейной связи между двумя
    переменными используется парная OLS-регрессия с HAC
    (Newey-West) стандартными ошибками.

    Поправка Холма применяется к совокупности парных тестов
    внутри одной корреляционной матрицы.

    Parameters
    ----------
    df : pd.DataFrame
        Таблица логарифмических доходностей.
    hac_maxlags : int
        Максимальное число лагов для HAC-оценки.

    Returns
    -------
    corr : pd.DataFrame
        Матрица корреляций Пирсона.
    p_raw : pd.DataFrame
        Необработанные HAC p-value.
    p_adj : pd.DataFrame
        HAC p-value с поправкой Холма.
    """

    if hac_maxlags < 0:
        raise ValueError("hac_maxlags не может быть отрицательным.")

    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.replace([np.inf, -np.inf], np.nan)

    columns = df.columns

    # Инициализация матриц
    corr = pd.DataFrame(
        np.nan,
        index=columns,
        columns=columns
    )

    p_raw = pd.DataFrame(
        np.nan,
        index=columns,
        columns=columns
    )

    p_adj = pd.DataFrame(
        np.nan,
        index=columns,
        columns=columns
    )

    # Диагональ корреляционной матрицы
    for col in columns:
        corr.loc[col, col] = 1.0

    valid_pairs = []

    # Перебор уникальных пар переменных
    for col_x, col_y in combinations(columns, 2):

        # Используем совместные наблюдения пары
        pair = df[[col_x, col_y]].dropna()

        x = pair[col_x]
        y = pair[col_y]

        # Проверка достаточности наблюдений
        if len(pair) <= hac_maxlags + 2:
            continue

        # Проверка наличия вариации в обоих рядах
        if x.nunique() < 2 or y.nunique() < 2:
            continue

        # Коэффициент корреляции Пирсона
        r = x.corr(y)

        corr.loc[col_x, col_y] = r
        corr.loc[col_y, col_x] = r

        # Парная линейная регрессия:
        # y = alpha + beta * x + error
        X = sm.add_constant(x.rename("x"))

        model = sm.OLS(y, X).fit(
            cov_type="HAC",
            cov_kwds={"maxlags": hac_maxlags},
            use_t=True
        )

        # H0: beta = 0
        p_value = model.pvalues["x"]

        p_raw.loc[col_x, col_y] = p_value
        p_raw.loc[col_y, col_x] = p_value

        valid_pairs.append(
            (col_x, col_y, p_value)
        )

    # Поправка Холма на множественные проверки
    if valid_pairs:

        raw_pvalues = [
            item[2] for item in valid_pairs
        ]

        adjusted_pvalues = multipletests(
            raw_pvalues,
            method="holm"
        )[1]

        for (col_x, col_y, _), p_value in zip(
            valid_pairs,
            adjusted_pvalues
        ):

            p_adj.loc[col_x, col_y] = p_value
            p_adj.loc[col_y, col_x] = p_value

    return corr, p_raw, p_adj

# Тепловые карты корреляционных матриц
FIGURES_DIR = BASE_DIR / "outputs" / "figures"
CORR_DIR = FIGURES_DIR / "correlation_matrices"

CORR_DIR.mkdir(parents=True, exist_ok=True)


def plot_corr_heatmap(
    corr: pd.DataFrame,
    title: str,
    save_path: Path
):
    """
    Строит тепловую карту корреляционной матрицы.

    Используется фиксированная цветовая шкала [-1; 1],
    обеспечивающая сопоставимость графиков.
    """

    fig, ax = plt.subplots(figsize=(10, 8))

    # Маскируем отсутствующие значения
    corr_values = np.ma.masked_invalid(corr.values)

    im = ax.imshow(
        corr_values,
        cmap="RdBu_r",
        vmin=-1,
        vmax=1
    )

    # Цветовая шкала
    cbar = fig.colorbar(
        im,
        ax=ax,
        ticks=np.linspace(-1, 1, 9)
    )

    cbar.set_label(
        "Pearson correlation coefficient",
        fontsize=10
    )

    # Подписи осей
    ax.set_xticks(range(len(corr.columns)))
    ax.set_yticks(range(len(corr.columns)))

    ax.set_xticklabels(
        corr.columns,
        rotation=45,
        ha="right"
    )

    ax.set_yticklabels(corr.columns)

    # Значения коэффициентов внутри ячеек
    for i in range(len(corr.columns)):
        for j in range(len(corr.columns)):

            value = corr.iloc[i, j]

            if pd.notna(value):

                text_color = (
                    "white" if abs(value) > 0.55 else "black"
                )

                ax.text(
                    j,
                    i,
                    f"{value:.2f}",
                    ha="center",
                    va="center",
                    color=text_color,
                    fontsize=8
                )

    ax.set_title(title, fontsize=12, pad=15)
    ax.set_aspect("equal")

    fig.tight_layout()

    # Сохранение изображения
    fig.savefig(
        save_path,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close(fig)

# Расчёт корреляционных матриц и построение графиков
CORRELATIONS = {}
P_VALUES = {}
ADJUSTED_P_VALUES = {}

for name, df in RETURNS.items():

    corr, p_raw, p_adj = corr_with_pvalues(
        df,
        hac_maxlags=5
    )

    # Сохраняем результаты в словарях
    CORRELATIONS[name] = corr
    P_VALUES[name] = p_raw
    ADJUSTED_P_VALUES[name] = p_adj

    print(
        f"\n========== Корреляционная матрица: {name} =========="
    )

    print("\nКоэффициенты корреляции Пирсона:")
    print(corr.round(3))

    print("\nHAC p-value:")
    print(p_raw.round(4))

    print("\nHAC p-value с поправкой Холма:")
    print(p_adj.round(4))

    # Сохраняем тепловую карту
    save_path = CORR_DIR / f"correlation_{name}.png"

    plot_corr_heatmap(
        corr,
        title=f"Correlation Matrix: {name}",
        save_path=save_path
    )


print(
    f"\nТепловые карты сохранены в папке: "
    f"{CORR_DIR.resolve()}"
)

# Выявление структурных сдвигов и деление на подпериоды
# Сохранение датасетов в csv для R
# Путь к папке для результатов теста Бай-Перрона
BP_DIR = BASE_DIR / 'outputs' / 'bai-perron'
BP_DIR.mkdir(parents=True, exist_ok=True)

EXPECTED_ASSETS = [
    'Gold',
    'Silver',
    'Cons_discr',
    'Cons_staples',
    'Energy',
    'Financials',
    'Inf_tech',
    'Health',
    'Industrials',
    'Materials',
    'Real_estate',
    'Telecom',
    'Utilities'
]

# ------------------------------------------------------------
# 2. Расчет логарифмических доходностей каждого актива
# ------------------------------------------------------------

returns_by_asset = {}

for asset in EXPECTED_ASSETS:

    # Извлекаем цены и приводим их к числовому формату
    prices = pd.to_numeric(
        series[asset][asset],
        errors='coerce'
    )

    # Сортировка дат по возрастанию
    prices = prices.sort_index()

    # Проверка уникальности дат
    if prices.index.has_duplicates:
        raise ValueError(
            f'Обнаружены повторяющиеся даты для актива {asset}'
        )

    # Удаление бесконечных и некорректных значений
    prices = prices.replace([np.inf, -np.inf], np.nan)

    # Цены должны быть положительными для расчета логарифма
    prices = prices.where(prices > 0)

    # Логарифмическая доходность, умноженная на 100
    returns = np.log(prices).diff() * 100

    returns.name = asset

    returns_by_asset[asset] = returns

# ------------------------------------------------------------
# 3. Объединение доходностей всех 13 активов
# ------------------------------------------------------------

# ВАЖНО:
# Используется outer join по датам.
# Индивидуальные торговые календари активов сохраняются.
# Пропуски не удаляются на уровне всей таблицы.

final_returns_all = pd.concat(
    returns_by_asset.values(),
    axis=1,
    join='outer'
).sort_index()

final_returns_all.index.name = 'Date'

# ------------------------------------------------------------
# 4. Контроль корректности итогового датасета
# ------------------------------------------------------------

# Проверка количества активов
assert final_returns_all.shape[1] == 13, (
    'Количество столбцов не равно 13'
)

# Проверка уникальности названий столбцов
assert final_returns_all.columns.is_unique, (
    'В датасете обнаружены повторяющиеся названия столбцов'
)

# Проверка наличия всех ожидаемых активов
assert set(final_returns_all.columns) == set(EXPECTED_ASSETS)

# Проверка хронологического порядка дат
assert final_returns_all.index.is_monotonic_increasing

# Проверка отсутствия бесконечных значений
assert not np.isinf(final_returns_all.select_dtypes(
    include=[np.number]
).to_numpy()).any()

# ------------------------------------------------------------
# 5. Сохранение CSV
# ------------------------------------------------------------

output_file = BP_DIR / 'final_returns_all.csv'

final_returns_all.to_csv(
    output_file,
    index_label='Date',
    date_format='%Y-%m-%d',
    encoding='utf-8'
)

# Разбиение на подпериоды
period1_c = final_returns_c.loc['2012-01-05':'2015-07-31']
period2_c = final_returns_c.loc['2015-08-01':'2016-12-31']
period3_c = final_returns_c.loc['2017-01-01':'2019-12-31']
period4_c = final_returns_c.loc['2020-01-01':'2020-12-31']
period5_c = final_returns_c.loc['2021-01-01':'2022-02-22']
period6_c = final_returns_c.loc['2022-02-23':'2022-12-31']
period7_c = final_returns_c.loc['2023-01-01':'2026-05-26']

levels1_c = final_levels_c.loc['2012-01-05':'2015-07-31']
levels2_c = final_levels_c.loc['2015-08-01':'2016-12-31']
levels3_c = final_levels_c.loc['2017-01-01':'2019-12-31']
levels4_c = final_levels_c.loc['2020-01-01':'2020-12-31']
levels5_c = final_levels_c.loc['2021-01-01':'2022-02-22']
levels6_c = final_levels_c.loc['2022-02-23':'2022-12-31']
levels7_c = final_levels_c.loc['2023-01-01':'2026-05-26']

period1_d = final_returns_d.loc['2012-01-05':'2015-07-31']
period2_d = final_returns_d.loc['2015-08-01':'2016-12-31']
period3_d = final_returns_d.loc['2017-01-01':'2019-12-31']
period4_d = final_returns_d.loc['2020-01-01':'2020-12-31']
period5_d = final_returns_d.loc['2021-01-01':'2022-02-22']
period6_d = final_returns_d.loc['2022-02-23':'2022-12-31']
period7_d = final_returns_d.loc['2023-01-01':'2026-05-26']

levels1_d = final_levels_d.loc['2012-01-05':'2015-07-31']
levels2_d = final_levels_d.loc['2015-08-01':'2016-12-31']
levels3_d = final_levels_d.loc['2017-01-01':'2019-12-31']
levels4_d = final_levels_d.loc['2020-01-01':'2020-12-31']
levels5_d = final_levels_d.loc['2021-01-01':'2022-02-22']
levels6_d = final_levels_d.loc['2022-02-23':'2022-12-31']
levels7_d = final_levels_d.loc['2023-01-01':'2026-05-26']

period1_m = final_returns_m.loc['2012-01-05':'2015-07-31']
period2_m = final_returns_m.loc['2015-08-01':'2016-12-31']
period3_m = final_returns_m.loc['2017-01-01':'2019-12-31']
period4_m = final_returns_m.loc['2020-01-01':'2020-12-31']
period5_m = final_returns_m.loc['2021-01-01':'2022-02-22']
period6_m = final_returns_m.loc['2022-02-23':'2022-12-31']
period7_m = final_returns_m.loc['2023-01-01':'2026-05-26']

levels1_m = final_levels_m.loc['2012-01-05':'2015-07-31']
levels2_m = final_levels_m.loc['2015-08-01':'2016-12-31']
levels3_m = final_levels_m.loc['2017-01-01':'2019-12-31']
levels4_m = final_levels_m.loc['2020-01-01':'2020-12-31']
levels5_m = final_levels_m.loc['2021-01-01':'2022-02-22']
levels6_m = final_levels_m.loc['2022-02-23':'2022-12-31']
levels7_m = final_levels_m.loc['2023-01-01':'2026-05-26']

# ============================================================
# Сохранение всех датасетов по подпериодам в CSV для R
# ============================================================
# Корневая папка проекта Sector Stocks
BASE_DIR = Path(__file__).resolve().parent

# Папка для подготовленных данных
DATA_FOR_MODELS_DIR = BASE_DIR / 'outputs' / 'data_for_models'

# Создание папки, если она отсутствует
DATA_FOR_MODELS_DIR.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# 1. Формирование словаря всех датасетов
# ------------------------------------------------------------

DATASETS_FOR_MODELS = {}

# Логарифмические доходности
for group in ['c', 'd', 'm']:
    for period in range(1, 8):
        dataset_name = f'period{period}_{group}'
        DATASETS_FOR_MODELS[dataset_name] = globals()[dataset_name]

# Логарифмические уровни
for group in ['c', 'd', 'm']:
    for period in range(1, 8):
        dataset_name = f'levels{period}_{group}'
        DATASETS_FOR_MODELS[dataset_name] = globals()[dataset_name]


# ------------------------------------------------------------
# 2. Проверка и сохранение датасетов
# ------------------------------------------------------------

for dataset_name, df in DATASETS_FOR_MODELS.items():

    # Проверка, что объект является DataFrame
    if not isinstance(df, pd.DataFrame):
        raise TypeError(
            f'{dataset_name} не является pandas DataFrame'
        )

    # Проверка уникальности дат
    if df.index.has_duplicates:
        raise ValueError(
            f'В датасете {dataset_name} обнаружены повторяющиеся даты'
        )

    # Проверка хронологического порядка
    if not df.index.is_monotonic_increasing:
        raise ValueError(
            f'Даты в {dataset_name} расположены не по возрастанию'
        )

    # Проверка отсутствия пустых значений
    if df.isna().any().any():
        raise ValueError(
            f'В датасете {dataset_name} обнаружены пропущенные значения'
        )

    # Сохранение в CSV
    output_file = DATA_FOR_MODELS_DIR / f'{dataset_name}.csv'

    df.to_csv(
        output_file,
        index_label='Date',
        date_format='%Y-%m-%d',
        encoding='utf-8'
    )

    print(
        f'{dataset_name}: сохранён, '
        f'{df.shape[0]} наблюдений, '
        f'{df.shape[1]} переменных'
    )


# ------------------------------------------------------------
# 3. Итоговая проверка экспорта
# ------------------------------------------------------------

saved_files = list(DATA_FOR_MODELS_DIR.glob('*.csv'))

print('\n======================================')
print('ЭКСПОРТ ЗАВЕРШЁН')
print('======================================')

print(f'Количество датасетов: {len(saved_files)}')
print(f'Папка сохранения: {DATA_FOR_MODELS_DIR.resolve()}')

if len(saved_files) != 42:
    raise RuntimeError(
        f'Ожидалось 42 CSV-файла, сохранено {len(saved_files)}'
    )

print('Все 42 датасета успешно сохранены.')
