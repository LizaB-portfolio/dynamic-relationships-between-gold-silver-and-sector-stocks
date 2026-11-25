
# ============================================================
# BAI-PERRON MULTIPLE STRUCTURAL BREAK ANALYSIS
# Sector Stocks Project
# ============================================================


# ------------------------------------------------------------
# 1. Установка и загрузка библиотек
# ------------------------------------------------------------

required_packages <- c(
  "here",
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "zoo",
  "strucchange",
  "purrr",
  "tibble"
)

# Установка отсутствующих пакетов
new_packages <- required_packages[
  !required_packages %in% installed.packages()[, "Package"]
]

if (length(new_packages) > 0) {
  install.packages(new_packages)
}

# Загрузка библиотек
invisible(
  lapply(required_packages, library, character.only = TRUE)
)


# ------------------------------------------------------------
# 2. Пути к данным и выходным папкам
# ------------------------------------------------------------

# Корневая директория проекта
PROJECT_ROOT <- here::here()

# Папка с исходным CSV
BP_DIR <- file.path(
  PROJECT_ROOT,
  "outputs",
  "bai-perron"
)

INPUT_FILE <- file.path(
  BP_DIR,
  "final_returns_all.csv"
)

# Папки для результатов
RESULTS_DIR <- file.path(
  BP_DIR,
  "results"
)

TABLES_DIR <- file.path(
  RESULTS_DIR,
  "tables"
)

FIGURES_DIR <- file.path(
  RESULTS_DIR,
  "figures"
)

MEAN_FIGURES_DIR <- file.path(
  FIGURES_DIR,
  "mean_breaks"
)

VARIANCE_FIGURES_DIR <- file.path(
  FIGURES_DIR,
  "variance_breaks"
)

BIC_FIGURES_DIR <- file.path(
  FIGURES_DIR,
  "bic"
)

# Создание директорий
output_dirs <- c(
  RESULTS_DIR,
  TABLES_DIR,
  FIGURES_DIR,
  MEAN_FIGURES_DIR,
  VARIANCE_FIGURES_DIR,
  BIC_FIGURES_DIR
)

for (dir in output_dirs) {
  dir.create(
    dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ------------------------------------------------------------
# 3. Параметры теста Бай — Перрона
# ------------------------------------------------------------

# Минимальная длина сегмента:
# 5% от общего количества наблюдений
MIN_SEGMENT <- 0.05

# Максимальное количество структурных сдвигов,
# рассматриваемых алгоритмом
MAX_BREAKS <- 8L

# Окно для расчета скользящего стандартного отклонения
ROLLING_WINDOW <- 21L


# ------------------------------------------------------------
# 4. Загрузка данных
# ------------------------------------------------------------

if (!file.exists(INPUT_FILE)) {
  stop(
    paste0(
      "Файл не найден: ",
      INPUT_FILE,
      "\nСначала запусти Python-код формирования CSV."
    )
  )
}

dataBP <- readr::read_csv(
  INPUT_FILE,
  show_col_types = FALSE,
  na = c("", "NA", "NaN")
)

# Проверка наличия даты
if (!"Date" %in% names(dataBP)) {
  stop("В датасете отсутствует столбец Date.")
}

# Преобразование Date в формат даты
dataBP <- dataBP %>%
  mutate(Date = as.Date(Date)) %>%
  arrange(Date)

# Проверка уникальности дат
if (anyDuplicated(dataBP$Date) > 0) {
  stop("В датасете обнаружены повторяющиеся даты.")
}

# Проверка количества активов
EXPECTED_ASSETS <- c(
  "Gold",
  "Silver",
  "Cons_discr",
  "Cons_staples",
  "Energy",
  "Financials",
  "Inf_tech",
  "Health",
  "Industrials",
  "Materials",
  "Real_estate",
  "Telecom",
  "Utilities"
)

missing_assets <- setdiff(
  EXPECTED_ASSETS,
  names(dataBP)
)

if (length(missing_assets) > 0) {
  stop(
    paste(
      "В датасете отсутствуют активы:",
      paste(missing_assets, collapse = ", ")
    )
  )
}

# Оставляем только необходимые столбцы
dataBP <- dataBP %>%
  select(Date, all_of(EXPECTED_ASSETS))

# Преобразование значений в числовой формат
dataBP <- dataBP %>%
  mutate(
    across(
      all_of(EXPECTED_ASSETS),
      as.numeric
    )
  )

cat("\n===== DATASET INFORMATION =====\n")

cat("Количество календарных дат:", nrow(dataBP), "\n")
cat("Количество активов:", length(EXPECTED_ASSETS), "\n")
cat("Первая дата:", as.character(min(dataBP$Date)), "\n")
cat("Последняя дата:", as.character(max(dataBP$Date)), "\n")


# ------------------------------------------------------------
# 5. Функция проведения теста Бай — Перрона
# ------------------------------------------------------------

run_bai_perron <- function(
    data,
    asset,
    model_type,
    h = MIN_SEGMENT,
    max_breaks = MAX_BREAKS
) {
  
  # Подготовка индивидуального ряда
  asset_data <- data %>%
    select(
      Date,
      Return = all_of(asset)
    ) %>%
    filter(
      !is.na(Date),
      is.finite(Return)
    ) %>%
    arrange(Date)
  
  # Проверка количества наблюдений
  if (nrow(asset_data) < 50) {
    warning(
      paste(
        "Недостаточно наблюдений для актива:",
        asset
      )
    )
    
    return(NULL)
  }
  
  # Формула модели
  if (model_type == "mean") {
    
    # Структурные сдвиги среднего
    formula <- Return ~ 1
    
  } else if (model_type == "variance") {
    
    # Структурные сдвиги второго момента доходности
    formula <- I(Return^2) ~ 1
    
  } else {
    stop("Неизвестный тип модели.")
  }
  
  # ----------------------------------------------------------
  # 5.1. Оценивание модели Бай — Перрона
  # ----------------------------------------------------------
  
  bp_model <- strucchange::breakpoints(
    formula,
    data = asset_data,
    h = h,
    breaks = max_breaks
  )
  
  # ----------------------------------------------------------
  # 5.2. Выбор количества структурных сдвигов по BIC
  # ----------------------------------------------------------
  
  bic_values <- AIC(
    bp_model,
    k = log(bp_model$nobs)
  )
  
  bic_values <- as.numeric(bic_values)
  
  bic_table <- tibble(
    Asset = asset,
    Model = model_type,
    Number_of_breaks = seq_along(bic_values) - 1L,
    BIC = bic_values
  )
  
  # Оптимальное количество структурных сдвигов
  selected_n <- which.min(bic_values) - 1L
  
  # ----------------------------------------------------------
  # 5.3. Извлечение дат структурных сдвигов
  # ----------------------------------------------------------
  
  if (selected_n == 0L) {
    
    break_idx <- integer(0)
    
  } else {
    
    selected_model <- strucchange::breakpoints(
      bp_model,
      breaks = selected_n
    )
    
    break_idx <- as.integer(
      selected_model$breakpoints
    )
    
    break_idx <- break_idx[
      !is.na(break_idx)
    ]
  }
  
  # ----------------------------------------------------------
  # 5.4. Формирование таблицы дат структурных сдвигов
  # ----------------------------------------------------------
  
  if (length(break_idx) > 0) {
    
    # break_idx — последнее наблюдение предыдущего сегмента.
    # Следующее наблюдение — начало нового сегмента.
    
    break_dates <- tibble(
      Asset = asset,
      Model = model_type,
      Break_Number = seq_along(break_idx),
      Break_Index = break_idx,
      Date_Before_Break = asset_data$Date[break_idx],
      Date_After_Break = asset_data$Date[break_idx + 1L],
      Selected_Number_of_Breaks = selected_n,
      BIC = min(bic_values)
    )
    
  } else {
    
    break_dates <- tibble(
      Asset = character(),
      Model = character(),
      Break_Number = integer(),
      Break_Index = integer(),
      Date_Before_Break = as.Date(character()),
      Date_After_Break = as.Date(character()),
      Selected_Number_of_Breaks = integer(),
      BIC = double()
    )
  }
  
  # ----------------------------------------------------------
  # 5.5. Статистика по сегментам
  # ----------------------------------------------------------
  
  segment_points <- c(
    0L,
    break_idx,
    nrow(asset_data)
  )
  
  segment_summary <- purrr::map_dfr(
    seq_len(length(segment_points) - 1L),
    function(i) {
      
      start_idx <- segment_points[i] + 1L
      end_idx <- segment_points[i + 1L]
      
      segment_data <- asset_data[
        start_idx:end_idx,
        ,
        drop = FALSE
      ]
      
      tibble(
        Asset = asset,
        Model = model_type,
        Segment = i,
        Start_Date = min(segment_data$Date),
        End_Date = max(segment_data$Date),
        Observations = nrow(segment_data),
        Mean_Return = mean(segment_data$Return),
        SD_Return = sd(segment_data$Return),
        Mean_Squared_Return = mean(
          segment_data$Return^2
        )
      )
    }
  )
  
  # ----------------------------------------------------------
  # 5.6. Подготовка данных для графиков
  # ----------------------------------------------------------
  
  if (model_type == "mean") {
    
    plot_data <- asset_data
    
    y_label <- "Daily log return, %"
    
    plot_title <- paste(
      "Bai-Perron structural breaks in mean:",
      asset
    )
    
    plot_data$Plot_Value <- plot_data$Return
    
  } else {
    
    # Скользящее стандартное отклонение
    plot_data <- asset_data %>%
      mutate(
        Plot_Value = zoo::rollapplyr(
          Return,
          width = ROLLING_WINDOW,
          FUN = sd,
          fill = NA,
          partial = FALSE
        )
      )
    
    y_label <- "Rolling standard deviation"
    
    plot_title <- paste(
      "Bai-Perron structural breaks in volatility:",
      asset
    )
  }
  
  # График временного ряда
  time_plot <- ggplot(
    plot_data,
    aes(
      x = Date,
      y = Plot_Value
    )
  ) +
    geom_line(
      linewidth = 0.35,
      na.rm = TRUE
    ) +
    labs(
      title = plot_title,
      x = "Date",
      y = y_label
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(
        size = 11,
        face = "bold"
      )
    )
  
  # Добавление вертикальных линий на датах начала новых режимов
  if (length(break_idx) > 0) {
    
    break_date_values <- asset_data$Date[
      break_idx + 1L
    ]
    
    time_plot <- time_plot +
      geom_vline(
        xintercept = as.numeric(break_date_values),
        linetype = "dashed",
        linewidth = 0.6
      )
  }
  
  # ----------------------------------------------------------
  # 5.7. График BIC
  # ----------------------------------------------------------
  
  bic_plot <- ggplot(
    bic_table,
    aes(
      x = Number_of_breaks,
      y = BIC
    )
  ) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2) +
    geom_vline(
      xintercept = selected_n,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    scale_x_continuous(
      breaks = bic_table$Number_of_breaks
    ) +
    labs(
      title = paste(
        "BIC selection:",
        asset,
        "-",
        model_type
      ),
      x = "Number of structural breaks",
      y = "BIC"
    ) +
    theme_minimal()
  
  # ----------------------------------------------------------
  # 5.8. Сохранение графиков
  # ----------------------------------------------------------
  
  safe_asset_name <- make.names(asset)
  
  if (model_type == "mean") {
    
    plot_dir <- MEAN_FIGURES_DIR
    
  } else {
    
    plot_dir <- VARIANCE_FIGURES_DIR
  }
  
  ggsave(
    filename = file.path(
      plot_dir,
      paste0(safe_asset_name, "_", model_type, "_breaks.png")
    ),
    plot = time_plot,
    width = 12,
    height = 6,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(
      BIC_FIGURES_DIR,
      paste0(safe_asset_name, "_", model_type, "_BIC.png")
    ),
    plot = bic_plot,
    width = 9,
    height = 5,
    dpi = 300
  )
  
  # ----------------------------------------------------------
  # 5.9. Итоговая информация
  # ----------------------------------------------------------
  
  cat(
    "\nAsset:", asset,
    "\nModel:", model_type,
    "\nObservations:", nrow(asset_data),
    "\nSelected number of breaks:", selected_n,
    "\nMinimum BIC:", min(bic_values),
    "\n"
  )
  
  # Возвращение результатов
  return(
    list(
      break_dates = break_dates,
      bic = bic_table,
      segment_summary = segment_summary,
      selected_n = selected_n,
      break_idx = break_idx
    )
  )
}


# ------------------------------------------------------------
# 6. Проведение теста для всех 13 активов
# ------------------------------------------------------------

all_results <- list()

for (asset in EXPECTED_ASSETS) {
  
  # Модель структурных сдвигов среднего
  mean_result <- run_bai_perron(
    data = dataBP,
    asset = asset,
    model_type = "mean"
  )
  
  # Модель структурных сдвигов второго момента
  variance_result <- run_bai_perron(
    data = dataBP,
    asset = asset,
    model_type = "variance"
  )
  
  all_results[[asset]] <- list(
    mean = mean_result,
    variance = variance_result
  )
}


# ------------------------------------------------------------
# 7. Объединение результатов всех активов
# ------------------------------------------------------------

# Даты структурных сдвигов
all_break_dates <- purrr::map_dfr(
  all_results,
  function(x) {
    purrr::map_dfr(
      x,
      function(model_result) {
        if (is.null(model_result)) {
          return(tibble())
        }
        
        model_result$break_dates
      }
    )
  }
)

# Значения BIC для всех моделей
all_bic <- purrr::map_dfr(
  all_results,
  function(x) {
    purrr::map_dfr(
      x,
      function(model_result) {
        if (is.null(model_result)) {
          return(tibble())
        }
        
        model_result$bic
      }
    )
  }
)

# Статистика по сегментам
all_segment_summary <- purrr::map_dfr(
  all_results,
  function(x) {
    purrr::map_dfr(
      x,
      function(model_result) {
        if (is.null(model_result)) {
          return(tibble())
        }
        
        model_result$segment_summary
      }
    )
  }
)


# ------------------------------------------------------------
# 8. Сохранение итоговых таблиц
# ------------------------------------------------------------

# 8.1. Даты структурных сдвигов
readr::write_csv(
  all_break_dates,
  file.path(
    TABLES_DIR,
    "break_dates.csv"
  )
)

# 8.2. Значения BIC и выбранное количество сдвигов
readr::write_csv(
  all_bic,
  file.path(
    TABLES_DIR,
    "bic_selection.csv"
  )
)

# 8.3. Описательная статистика по сегментам
readr::write_csv(
  all_segment_summary,
  file.path(
    TABLES_DIR,
    "segment_summary.csv"
  )
)


# ------------------------------------------------------------
# 9. Итоговая сводка
# ------------------------------------------------------------

selected_breaks_summary <- all_break_dates %>%
  group_by(
    Asset,
    Model
  ) %>%
  summarise(
    Number_of_Breaks = n(),
    Break_Dates = paste(
      as.character(Date_After_Break),
      collapse = "; "
    ),
    .groups = "drop"
  )

cat("\n\n===== BAI-PERRON RESULTS =====\n")

print(selected_breaks_summary)

cat("\nРезультаты сохранены в:\n")
cat(RESULTS_DIR, "\n")

cat("\n===== ANALYSIS COMPLETED =====\n")

