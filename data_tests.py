
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import statsmodels

from statsmodels.tsa.stattools import adfuller, kpss
from statsmodels.tsa.vector_ar.var_model import VAR

# ============================================================
# 1. НАСТРОЙКИ
# ============================================================

# Корневая папка проекта определяется автоматически:
# скрипт должен находиться в корне Sector Stocks.
BASE_DIR = Path(__file__).resolve().parent

INPUT_DIR = BASE_DIR / "outputs" / "data_for_models"
OUTPUT_DIR = BASE_DIR / "outputs" / "econometric_tests"
TABLES_DIR = OUTPUT_DIR / "tables"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
TABLES_DIR.mkdir(parents=True, exist_ok=True)

ALPHA = 0.05

# ADF
ADF_MAXLAG = 12
ADF_AUTOLAG = "BIC"

# Выбор порядка VAR
VAR_MAXLAGS = 12
VAR_TREND = "c"

# Составы трёх групп моделей.
GROUPS = {
    "c": [
        "Gold",
        "Silver",
        "Cons_discr",
        "Financials",
        "Industrials",
        "Materials",
    ],
    "d": [
        "Gold",
        "Silver",
        "Cons_staples",
        "Health",
        "Utilities",
    ],
    "m": [
        "Gold",
        "Silver",
        "Inf_tech",
        "Real_estate",
        "Telecom",
        "Energy",
    ],
}

PERIODS = range(1, 8)
GROUP_CODES = tuple(GROUPS.keys())


# ============================================================
# 2. ФОРМИРОВАНИЕ ИМЁН ФАЙЛОВ
# ============================================================

def expected_dataset_names():
    """
    Формирует список 42 датасетов:
    21 датасет доходностей и 21 датасет логарифмических уровней.
    """
    names = []

    for prefix in ("period", "levels"):
        for group in GROUP_CODES:
            for period in PERIODS:
                names.append(f"{prefix}{period}_{group}")

    return names


def get_expected_columns(dataset_name):
    group = dataset_name.rsplit("_", 1)[1]
    return GROUPS[group]


# ============================================================
# 3. ЗАГРУЗКА И ПРОВЕРКА ДАННЫХ
# ============================================================

def load_datasets():
    datasets = {}

    if not INPUT_DIR.exists():
        raise FileNotFoundError(
            f"Не найдена папка с исходными CSV-файлами: {INPUT_DIR}"
        )

    for dataset_name in expected_dataset_names():

        file_path = INPUT_DIR / f"{dataset_name}.csv"

        if not file_path.exists():
            raise FileNotFoundError(
                f"Не найден файл: {file_path}"
            )

        df = pd.read_csv(
            file_path,
            parse_dates=["Date"],
            index_col="Date",
        )

        if df.empty:
            raise ValueError(
                f"Датасет {dataset_name} пуст."
            )

        if not isinstance(df.index, pd.DatetimeIndex):
            raise TypeError(
                f"В {dataset_name} индекс Date не распознан как дата."
            )

        if df.index.has_duplicates:
            raise ValueError(
                f"В {dataset_name} обнаружены повторяющиеся даты."
            )

        if not df.index.is_monotonic_increasing:
            raise ValueError(
                f"Даты в {dataset_name} не отсортированы."
            )

        expected = get_expected_columns(dataset_name)

        missing = sorted(set(expected) - set(df.columns))
        extra = sorted(set(df.columns) - set(expected))

        if missing or extra:
            raise ValueError(
                f"Неверный состав переменных в {dataset_name}. "
                f"Отсутствуют: {missing}; лишние: {extra}."
            )

        # Сохраняем установленный порядок переменных.
        df = df[expected].apply(pd.to_numeric, errors="coerce")

        if df.isna().any().any():
            raise ValueError(
                f"В {dataset_name} есть пропуски "
                f"или нечисловые значения."
            )

        if not np.isfinite(df.to_numpy(dtype=float)).all():
            raise ValueError(
                f"В {dataset_name} обнаружены inf или -inf."
            )

        # Проверка минимального размера выборки.
        if len(df) <= ADF_MAXLAG + 5:
            raise ValueError(
                f"Недостаточно наблюдений в {dataset_name}: "
                f"n={len(df)}."
            )

        datasets[dataset_name] = df

    if len(datasets) != 42:
        raise RuntimeError(
            f"Ожидалось 42 датасета, загружено {len(datasets)}."
        )

    return datasets


def get_metadata(dataset_name, df):
    prefix, group = dataset_name.split("_")

    period = int(
        prefix.replace("period", "").replace("levels", "")
    )

    return {
        "dataset": dataset_name,
        "data_type": (
            "returns" if prefix.startswith("period")
            else "log_levels"
        ),
        "period": period,
        "group": group,
        "n_obs": len(df),
        "n_variables": df.shape[1],
        "start_date": df.index.min().date().isoformat(),
        "end_date": df.index.max().date().isoformat(),
    }


# ============================================================
# 4. ADF-ТЕСТ
# ============================================================

def run_adf(series, regression):
    """
    H0: ряд имеет единичный корень.

    regression:
        'c'  — константа;
        'ct' — константа и линейный тренд.
    """

    result = adfuller(
        series.to_numpy(dtype=float),
        maxlag=ADF_MAXLAG,
        regression=regression,
        autolag=ADF_AUTOLAG,
    )

    statistic, pvalue, usedlag, nobs, critical_values, icbest = result # pyright: ignore[reportAssignmentType]

    return {
        "adf_statistic": statistic,
        "adf_pvalue": pvalue,
        "adf_used_lag": usedlag,
        "adf_nobs": nobs,
        "adf_critical_1pct": critical_values["1%"],
        "adf_critical_5pct": critical_values["5%"],
        "adf_critical_10pct": critical_values["10%"],
        "adf_icbest": icbest,
    }


# ============================================================
# 5. KPSS-ТЕСТ
# ============================================================

def run_kpss(series, regression):
    """
    H0: ряд стационарен относительно:
        'c'  — константы;
        'ct' — детерминированного линейного тренда.
    """

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")

        statistic, pvalue, used_lags, critical_values = kpss(
            series.to_numpy(dtype=float),
            regression=regression,
            nlags="auto",
        )

    warning_text = " | ".join(
        str(w.message) for w in caught
    )

    return {
        "kpss_statistic": statistic,
        "kpss_pvalue": pvalue,
        "kpss_used_lags": used_lags,
        "kpss_critical_1pct": critical_values["1%"],
        "kpss_critical_5pct": critical_values["5%"],
        "kpss_critical_10pct": critical_values["10%"],
        "kpss_warning": warning_text,
    }


# ============================================================
# 6. ADF И KPSS ДЛЯ 42 ДАТАСЕТОВ
# ============================================================

def run_stationarity_tests(datasets):
    rows = []

    for dataset_name, df in datasets.items():

        is_levels = dataset_name.startswith("levels")

        # Для логарифмов цен и индексов:
        # тестируем стационарность относительно тренда.
        #
        # Для доходностей:
        # тестируем стационарность относительно константы.
        regression = "ct" if is_levels else "c"

        for column in df.columns:

            series = df[column]

            adf_result = run_adf(
                series,
                regression=regression,
            )

            kpss_result = run_kpss(
                series,
                regression=regression,
            )

            adf_reject = adf_result["adf_pvalue"] < ALPHA
            kpss_reject = kpss_result["kpss_pvalue"] < ALPHA

            # ADF отвергает H0 единичного корня:
            # свидетельство в пользу стационарности.
            #
            # KPSS отвергает H0 стационарности:
            # свидетельство против стационарности.
            if adf_reject and not kpss_reject:
                conclusion = "stationary"
            elif not adf_reject and kpss_reject:
                conclusion = "non_stationary"
            elif adf_reject and kpss_reject:
                conclusion = "conflicting_results"
            else:
                conclusion = "inconclusive"

            row = {
                **get_metadata(dataset_name, df),
                "variable": column,
                "regression": regression,
                "alpha": ALPHA,
                **adf_result,
                **kpss_result,
                "adf_reject_unit_root": adf_reject,
                "kpss_reject_stationarity": kpss_reject,
                "combined_conclusion": conclusion,
            }

            rows.append(row)

    return pd.DataFrame(rows)


# ============================================================
# 7. ВЫБОР ОПТИМАЛЬНОГО ПОРЯДКА VAR
# ============================================================

def get_effective_maxlags(df):
    """
    Ограничивает максимальный порядок VAR с учётом размера выборки
    и числа эндогенных переменных.

    Это предотвращает попытку оценить слишком много параметров
    относительно доступного числа наблюдений.
    """
    n_obs = len(df)
    n_vars = df.shape[1]

    # Консервативное ограничение по числу наблюдений
    # на уравнение VAR.
    sample_limit = n_obs // (n_vars + 1) - 1

    effective_maxlags = min(
        VAR_MAXLAGS,
        sample_limit,
    )

    if effective_maxlags < 1:
        raise ValueError(
            f"Недостаточно наблюдений для выбора лагов VAR: "
            f"n={n_obs}, variables={n_vars}."
        )

    return effective_maxlags


def run_lag_selection(datasets):
    summary_rows = []
    detail_rows = []

    for dataset_name, df in datasets.items():

        effective_maxlags = get_effective_maxlags(df)

        model = VAR(df)

        selection = model.select_order(
            maxlags=effective_maxlags,
            trend=VAR_TREND,
        )

        metadata = get_metadata(dataset_name, df)

        selected_orders = selection.selected_orders

        summary_row = {
            **metadata,
            "maxlags_tested": effective_maxlags,
            "trend": VAR_TREND,
            "selected_lag_aic": selected_orders.get("aic"),
            "selected_lag_bic": selected_orders.get("bic"),
            "selected_lag_hqic": selected_orders.get("hqic"),
            "selected_lag_fpe": selected_orders.get("fpe"),
        }

        summary_rows.append(summary_row)

        # Сохраняем значения критериев для каждого
        # проверенного порядка лага.
        for criterion, values in selection.ics.items():

            for lag, value in enumerate(values):

                detail_rows.append({
                    **metadata,
                    "criterion": criterion,
                    "lag_order": lag,
                    "criterion_value": float(value),
                })

    summary_df = pd.DataFrame(summary_rows)
    details_df = pd.DataFrame(detail_rows)

    return summary_df, details_df

# ============================================================
# 8. ПОДГОТОВКА ЧИТАЕМЫХ ТАБЛИЦ И ИХ СОХРАНЕНИЕ В PNG
# ============================================================

def render_table_png(df, title, file_stem, rows_per_image=35):
    """
    Сохраняет DataFrame в PNG.
    Длинные таблицы автоматически разбиваются на несколько изображений.
    """
    if df.empty:
        return

    for part_number, start in enumerate(
        range(0, len(df), rows_per_image),
        start=1,
    ):
        part = df.iloc[start:start + rows_per_image].copy()

        height = max(3.0, 1.1 + 0.36 * (len(part) + 1))
        width = max(12.0, min(28.0, 1.7 * len(part.columns)))

        fig, ax = plt.subplots(figsize=(width, height))
        ax.axis("off")

        display = part.copy().fillna("")

        table = ax.table(
            cellText=display.astype(str).values,
            colLabels=display.columns,
            cellLoc="center",
            loc="center",
        )

        table.auto_set_font_size(False)
        table.set_fontsize(7)
        table.scale(1.0, 1.35)

        if len(df) <= rows_per_image:
            figure_title = title
        else:
            figure_title = f"{title} (часть {part_number})"

        ax.set_title(
            figure_title,
            fontsize=11,
            pad=16,
        )

        fig.tight_layout()

        suffix = (
            f"_part{part_number}"
            if len(df) > rows_per_image
            else ""
        )

        fig.savefig(
            TABLES_DIR / f"{file_stem}{suffix}.png",
            dpi=220,
            bbox_inches="tight",
        )

        plt.close(fig)


def save_readable_table(df, title, file_stem, rows_per_image=35):
    """
    Сохраняет читаемую таблицу одновременно в CSV и PNG.
    """
    df.to_csv(
        TABLES_DIR / f"{file_stem}.csv",
        index=False,
        encoding="utf-8-sig",
    )

    render_table_png(
        df,
        title,
        file_stem,
        rows_per_image=rows_per_image,
    )


# ------------------------------------------------------------
# 8.1. Читаемая таблица ADF и KPSS
# ------------------------------------------------------------

def build_stationarity_readable(stationarity_results):
    """
    Формирует компактную таблицу с результатами ADF и KPSS.
    """
    columns = [
        "data_type",
        "group",
        "period",
        "variable",
        "n_obs",
        "start_date",
        "end_date",
        "regression",
        "adf_statistic",
        "adf_pvalue",
        "adf_used_lag",
        "kpss_statistic",
        "kpss_pvalue",
        "kpss_used_lags",
        "combined_conclusion",
    ]

    out = stationarity_results[columns].copy()

    out = out.rename(columns={
        "data_type": "Тип данных",
        "group": "Группа",
        "period": "Подпериод",
        "variable": "Переменная",
        "n_obs": "Наблюдения",
        "start_date": "Начало выборки",
        "end_date": "Конец выборки",
        "regression": "Детерминированные компоненты",
        "adf_statistic": "ADF: статистика",
        "adf_pvalue": "ADF: p-value",
        "adf_used_lag": "ADF: лаг",
        "kpss_statistic": "KPSS: статистика",
        "kpss_pvalue": "KPSS: p-value",
        "kpss_used_lags": "KPSS: лаги",
        "combined_conclusion": "Общий вывод",
    })

    return out.round(4)


# ------------------------------------------------------------
# 8.2. Матрицы итоговых статусов стационарности
# ------------------------------------------------------------

def build_stationarity_status_matrix(
    stationarity_results,
    data_type,
    group,
):
    subset = stationarity_results[
        (stationarity_results["data_type"] == data_type)
        & (stationarity_results["group"] == group)
    ].copy()

    if subset.empty:
        return pd.DataFrame()

    pivot = subset.pivot(
        index="variable",
        columns="period",
        values="combined_conclusion",
    )

    pivot = pivot.reindex(GROUPS[group])

    pivot.columns = [
        f"Подпериод {col}"
        for col in pivot.columns
    ]

    pivot.index.name = "Переменная"

    return pivot.reset_index()


# ------------------------------------------------------------
# 8.3. Читаемая таблица выбора лагов VAR
# ------------------------------------------------------------

def build_lag_readable(lag_summary):
    columns = [
        "data_type",
        "group",
        "period",
        "n_obs",
        "n_variables",
        "maxlags_tested",
        "selected_lag_aic",
        "selected_lag_bic",
        "selected_lag_hqic",
        "selected_lag_fpe",
    ]

    out = lag_summary[columns].copy()

    out = out.rename(columns={
        "data_type": "Тип данных",
        "group": "Группа",
        "period": "Подпериод",
        "n_obs": "Наблюдения",
        "n_variables": "Число переменных",
        "maxlags_tested": "Максимальный лаг в проверке",
        "selected_lag_aic": "Лаг AIC",
        "selected_lag_bic": "Лаг BIC",
        "selected_lag_hqic": "Лаг HQIC",
        "selected_lag_fpe": "Лаг FPE",
    })

    return out

# ------------------------------------------------------------
# 8.4. Сохранение всех читаемых таблиц
# ------------------------------------------------------------

def save_all_tables(
    stationarity_results,
    lag_summary,
):
    """
    Формирует читаемые CSV и PNG для:
    - ADF и KPSS;
    - итоговых статусов стационарности;
    - выбора порядка VAR.
    """

    # ADF + KPSS: общая таблица.
    stationarity_readable = build_stationarity_readable(
        stationarity_results
    )

    save_readable_table(
        stationarity_readable,
        "ADF и KPSS: результаты тестов стационарности",
        "stationarity_adf_kpss_all",
        rows_per_image=35,
    )

    # Матрицы статусов и подробные результаты отдельно
    # для логарифмических уровней и доходностей.
    for data_type, type_label in (
        ("log_levels", "levels"),
        ("returns", "returns"),
    ):
        for group in GROUP_CODES:

            matrix = build_stationarity_status_matrix(
                stationarity_results,
                data_type,
                group,
            )

            save_readable_table(
                matrix,
                f"Статусы стационарности: {type_label}, группа {group}",
                f"stationarity_matrix_{type_label}_{group}",
            )

            mask = (
                (stationarity_results["data_type"] == data_type)
                & (stationarity_results["group"] == group)
            )

            detail = stationarity_readable.loc[mask].copy()

            save_readable_table(
                detail,
                f"Детальные результаты ADF и KPSS: {type_label}, группа {group}",
                f"stationarity_detail_{type_label}_{group}",
            )

    # Выбор лагов VAR: общая таблица.
    lag_readable = build_lag_readable(lag_summary)

    save_readable_table(
        lag_readable,
        "Выбор порядка VAR по информационным критериям",
        "lag_selection_readable",
    )

    # Выбор лагов VAR отдельно по группам.
    for group in GROUP_CODES:

        group_table = lag_readable[
            lag_readable["Группа"] == group
        ]

        save_readable_table(
            group_table,
            f"Выбор лагов VAR: группа {group}",
            f"lag_selection_group_{group}",
        )

# ============================================================
# 10. ОСНОВНОЙ ЗАПУСК И СОХРАНЕНИЕ РЕЗУЛЬТАТОВ
# ============================================================

def main():

    print("Версия statsmodels:", statsmodels.__version__)
    print("Загрузка датасетов...")

    datasets = load_datasets()

    print(f"Загружено датасетов: {len(datasets)}")

    # Инвентаризация файлов.
    inventory = pd.DataFrame([
        get_metadata(name, df)
        for name, df in datasets.items()
    ])

    inventory.to_csv(
        OUTPUT_DIR / "dataset_inventory.csv",
        index=False,
        encoding="utf-8-sig",
    )

    # ADF и KPSS.
    print("Выполнение ADF и KPSS...")

    stationarity_results = run_stationarity_tests(datasets)

    stationarity_results.to_csv(
        OUTPUT_DIR / "stationarity_tests.csv",
        index=False,
        encoding="utf-8-sig",
    )

    # Выбор порядка VAR.
    print("Выбор порядка VAR по AIC, BIC, HQIC и FPE...")

    lag_summary, lag_details = run_lag_selection(datasets)

    lag_summary.to_csv(
        OUTPUT_DIR / "lag_selection_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )

    lag_details.to_csv(
        OUTPUT_DIR / "lag_selection_details.csv",
        index=False,
        encoding="utf-8-sig",
    )

    # Читаемые CSV и PNG-таблицы.
    print("Формирование читаемых таблиц и изображений...")

    save_all_tables(
        stationarity_results,
        lag_summary,
    )

    # Контроль ожидаемого количества результатов.
    expected_stationarity_rows = sum(
        len(GROUPS[group])
        for group in GROUP_CODES
    ) * len(PERIODS) * 2

    print("\nРабота завершена.")
    print("Папка результатов:", OUTPUT_DIR)
    print("Читаемые таблицы и изображения:", TABLES_DIR)

    print(
        "Строк в stationarity_tests:",
        len(stationarity_results),
        "| ожидалось:",
        expected_stationarity_rows,
    )

    print(
        "Датасетов в lag_selection_summary:",
        len(lag_summary),
        "| ожидалось: 42",
    )


if __name__ == "__main__":
    main()
