#' Convert a list of p-values to a long table
#'
#' @param pvalues_list A list recording the p-values from the SVG test.
#' Each element is a data.frame corresponding to a cell type.
#' @param which_pvalue A character specifying the column name of p-values. The default is `combined_pvalue`.
#'
#' @return A data.frame for the p-values in a long table.
#'
#' @export
mcubePvalueList2LongTable <- function(pvalues_list, which_pvalue = "combined_pvalue") {
  pvalues_df <- do.call(
    rbind,
    lapply(
      names(pvalues_list),
      function(x) {
        data.frame(
          celltype = x,
          gene = rownames(pvalues_list[[x]]),
          pvalue = pvalues_list[[x]][, which_pvalue]
        )
      }
    )
  )
  return(pvalues_df)
}

#' Generate a data.frame for QQ plot
#'
#' @param pvalues_df A data.frame with p-values.
#' @param which_pvalue A character string specifying the column name of p-values. The default is `combined_pvalue`.
#' @param minus_log10p_max A numeric value specifying the maximum value of -log10(p-value). The default is `5`.
#'
#' @return A data.frame for QQ plot.
#'
#' @export
mcubeQQPlotDF <- function(pvalues_df, which_pvalue = "combined_pvalue", minus_log10p_max = 5) {
  pvalues_order <- order(pvalues_df[, which_pvalue])
  pvalues <- pvalues_df[pvalues_order, which_pvalue]
  n_pvalues <- length(pvalues)
  pvalues_theoretical <- ppoints(n_pvalues)
  minus_log10p <- -log10(pvalues)
  minus_log10p[minus_log10p > minus_log10p_max] <- minus_log10p_max
  minus_log10p_theoretical <- -log10(pvalues_theoretical)

  qqplot_df <- data.frame(
    gene = rownames(pvalues_df)[pvalues_order],
    pvalue = pvalues,
    pvalue_theoretical = pvalues_theoretical,
    minus_log10p = minus_log10p,
    minus_log10p_theoretical = minus_log10p_theoretical
  )
  return(qqplot_df)
}

#' Rotate a 2D spatial coordinate matrix
#'
#' @param theta The angle of rotation in degrees.
#'
#' @return A rotation matrix of the given angle.
#'
#' @export
rotation_matrix_2d <- function(theta) {
  theta_rad <- theta * (pi / 180) # Convert degrees to radians
  matrix(
    c(
      cos(theta_rad), -sin(theta_rad),
      sin(theta_rad), cos(theta_rad)
    ),
    nrow = 2,
    ncol = 2,
    byrow = TRUE
  )
}

#' Format the GWAS summary statistics data
#'
#' @details
#' Format the GWAS summary statistics data for linkage disequilibrium score regression (LDSC) (\url{https://doi.org/10.1038/ng.3211}).
#' Credit goes to `ldsc` (\url{https://github.com/bulik/ldsc}) and the R package `MR-APSS` (\url{https://github.com/YangLabHKUST/MR-APSS}).
#'
#' @importFrom readr read_delim
#' @importFrom dplyr inner_join select filter distinct one_of
#' @importFrom tidyr drop_na
#' @importFrom stats qchisq pchisq sd median
#'
#' @param sumstats Summary statistics data or name of the summary statistics data file.
#' @param delim A single character used to separate fields within a record.
#' @param snps_merge A data.frame with SNPs to extract.
#' @param snps_mhc A set of SNPs that needed to be removed. For example, SNPs in MHC region.
#' @param dbSNP A tibble with two columns `rsID` = (rs) identifiers and `SNP` = chr:location. The data is from the dbSNP153 database (https://genome.ucsc.edu/FAQ/FAQdownloads.html#snp).
#' @param snp_col Name of a column with rs numbers. The default is `SNP`.
#' @param beta_col Name of a column with effect sizes (= log(or)). The default is `BETA`.
#' @param or_col Name of a column with odds ratios. The default is `OR`.
#' @param se_col Name of a column with standard errors. The default is `SE`.
#' @param freq_col Name of a column with effect allele frequencies. The default is `FREQ`.
#' @param a1_col Name of a column with effect alleles. The default is `A1`.
#' @param a2_col Name of a column with the other alleles. The default is `A2`.
#' @param p_col Name of a column with p-values. The default is `P`.
#' @param n_col Name of a column with sample sizes. The default is `N`.
#' @param ncase_col Name of a column with the number of cases. The default is `N_CASE`.
#' @param ncontrol_col Name of a column with the number of controls. The default is `N_CONTROL`.
#' @param n Sample size. If the column for sample size is not available, users can use argument "n" to specify the total sample size for each SNP.
#' @param z_col Name of a column with z scores. The default is `Z`.
#' @param info_col Name of a column with imputation INFO. The default is `INFO`.
#' @param log_pval Logical, whether the p-value is in -log10 scale. The default is `FALSE`.
#' @param chi2_max SNPs with tested chi^2 statistics large than chi2_max will be removed. The default is `80`.
#' @param min_freq SNPs with allele frequency less than min_freq will be removed. The default is `0.05`.
#'
#' @return Formatted summary statistics data
#'
#' @export
format_sumstats <- function(
    sumstats, delim = NULL,
    snps_merge = NULL, snps_mhc = NULL, dbSNP = NULL,
    snp_col = "SNP", beta_col = "BETA", or_col = "OR", se_col = "SE",
    freq_col = "FREQ", a1_col = "A1", a2_col = "A2", p_col = "P",
    n_col = "N", ncase_col = "N_CASE", ncontrol_col = "N_CONTROL", n = NULL,
    z_col = "Z", info_col = "INFO", log_pval = FALSE,
    chi2_max = NULL, min_freq = 0.05) {
  if (is.null(sumstats)) {
    stop("please provide the information on summary statistics!")
  }

  if (is.character(sumstats)) {
    message("Begin reading in summary statistics data...")
    sumstats <- readr::read_delim(sumstats, delim = delim, comment = "##")
  }

  message("Begin formatting...")
  message("The raw dataset has ", nrow(sumstats), " lines.")

  cols <- c(
    snp_col, beta_col, or_col, se_col,
    freq_col, a1_col, a2_col, p_col,
    ncase_col, ncontrol_col, n_col,
    z_col, info_col
  )
  cols <- names(sumstats)[names(sumstats) %in% cols]
  out <- c("SNP", "A1", "A2", "Z", "N", "CHI2", "P")

  sumstats <- sumstats[, cols]

  # check SNP
  if (!snp_col %in% names(sumstats)) {
    stop("SNP not found!")
  }

  # check A1/A2
  if ((!a1_col %in% names(sumstats)) | (!a2_col %in% names(sumstats))) {
    stop("A1/A2 not found!")
  }

  # check signed statistics
  if ((!beta_col %in% names(sumstats)) &
    (!or_col %in% names(sumstats)) &
    (!z_col %in% names(sumstats))) {
    stop("Signed statistics not found!")
  }

  # check sample size n
  if ((!n_col %in% names(sumstats)) & is.null(n) &
    (!(ncase_col %in% names(sumstats) & ncontrol_col %in% names(sumstats)))) {
    stop("information for sample size not found!")
  }

  # check sample size n
  if (ncase_col %in% names(sumstats)) {
    names(sumstats)[which(names(sumstats) == ncase_col)[1]] <- "N_CASE"
    if (!is.numeric(sumstats$N_CASE)) {
      message("n_case column is not numeric data. Coercing...")
      sumstats$N_CASE <- as.numeric(as.character(sumstats$N_CASE))
    }
  }

  if (ncontrol_col %in% names(sumstats)) {
    names(sumstats)[which(names(sumstats) == ncontrol_col)[1]] <- "N_CONTROL"
    if (!is.numeric(sumstats$N_CONTROL)) {
      message("n_control column is not numeric data. Coercing...")
      sumstats$N_CONTROL <- as.numeric(as.character(sumstats$N_CONTROL))
    }
  }

  if (n_col %in% names(sumstats)) {
    names(sumstats)[which(names(sumstats) == n_col)[1]] <- "N"
    if (!is.numeric(sumstats$N)) {
      message("Sample size column is not numeric data. Coercing...")
      sumstats$N <- as.numeric(as.character(sumstats$N))
    }
  } else if ("N_CASE" %in% names(sumstats) & "N_CONTROL" %in% names(sumstats)) {
    message("Generate sample size from sample size of case and control.")
    sumstats$N <- sumstats$N_CASE + sumstats$N_CONTROL
  } else if (!is.null(n)) {
    message("Generate sample size from specified sample size.")
    sumstats$N <- n
  }
  sumstats$N <- as.integer(sumstats$N)

  names(sumstats)[which(names(sumstats) == snp_col)[1]] <- "SNP"

  # Format SNP id
  sumstats$SNP <- tolower(sumstats$SNP)
  sumstats$SNP <- gsub("[[:space:]]", "", sumstats$SNP)
  if (grepl("^rs[0-9]+:", sumstats$SNP[1])) {
    sumstats$SNP <- gsub(":.*", "", sumstats$SNP)
  }

  # Merge with dbSNP
  if (!is.null(dbSNP)) {
    sumstats <- dplyr::inner_join(sumstats, dbSNP, by = "SNP")
    sumstats <- dplyr::select(sumstats, -("SNP"))
    names(sumstats)[which(names(sumstats) == "rsID")[1]] <- "SNP"
    message(
      "Merge the data with dbSNP..., remaining ",
      nrow(sumstats), " SNPs."
    )
  }

  # Remove NAs
  sumstats <- tidyr::drop_na(sumstats, SNP)
  message("Remove missing data..., remaining ", nrow(sumstats), " SNPs.")

  # Merge with HapMap 3 SNPs lists
  if (!is.null(snps_merge)) {
    snps_merge <- dplyr::select(snps_merge, SNP, A1, A2)
    names(snps_merge) <- c("SNP", "ref_A1", "ref_A2")
    sumstats <- dplyr::inner_join(sumstats, snps_merge, by = "SNP")
    message(
      "Merge the data with the HapMap 3 SNPs list by SNP..., remaining ",
      nrow(sumstats), " SNPs."
    )
  }

  # remove MHC SNPs
  if (!is.null(snps_mhc)) {
    sumstats <- dplyr::filter(sumstats, !SNP %in% snps_mhc)
    message("Remove the SNPs in MHC region..., remaining ", nrow(sumstats), " SNPs.")
  }

  # check imputation INFO
  if (info_col %in% names(sumstats)) {
    names(sumstats)[which(names(sumstats) == info_col)[1]] <- "INFO"
    sumstats <- dplyr::filter(sumstats, INFO > 0.9 | is.na(INFO))
    message(
      "Remove the SNPs with imputation INFO less than 0.9 (keeping NAs)..., remaining ",
      nrow(sumstats), " SNPs."
    )
  }

  # check effect allele (A1)
  if (a1_col %in% names(sumstats)) {
    names(sumstats)[which(names(sumstats) == a1_col)[1]] <- "A1"

    if (is.logical(sumstats$A1)) {
      sumstats$A1 <- substr(as.character(sumstats$A1), 1, 1)
    }

    if (!is.character(sumstats$A1)) {
      message("Effect allele column is not character data. Coercing...")
      sumstats$A1 <- as.character(sumstats$A1)
    }

    sumstats$A1 <- toupper(sumstats$A1)

    sumstats <- dplyr::filter(sumstats, grepl("^[ACTG]+$", A1))
    message(
      "Remove the SNPs whose effect allele has the value that is not A/C/T/G...., remaining ",
      nrow(sumstats), " SNPs."
    )
  }

  # return(sumstats)

  # check the other allele (A2)
  if (a2_col %in% names(sumstats)) {
    names(sumstats)[which(names(sumstats) == a2_col)[1]] <- "A2"

    if (is.logical(sumstats$A2)) {
      sumstats$A2 <- substr(as.character(sumstats$A2), 1, 1)
    }

    if (!is.character(sumstats$A2)) {
      message("The other allele column is not character data. Coercing...")
      sumstats$A2 <- as.character(sumstats$A2)
    }

    sumstats$A2 <- toupper(sumstats$A2)

    sumstats <- dplyr::filter(sumstats, grepl("^[ACTG]+$", A2))
    message(
      "Remove the SNPs whose non effect allele has the value that is not A/C/T/G...., remaining ",
      nrow(sumstats), " SNPs."
    )
  }

  # remove ambiguous SNPs
  # A-T, C-G
  sumstats <- dplyr::filter(
    sumstats,
    (A1 == "A" & A2 == "C") |
      (A1 == "A" & A2 == "G") |
      (A1 == "T" & A2 == "C") |
      (A1 == "T" & A2 == "G") |
      (A1 == "C" & A2 == "A") |
      (A1 == "C" & A2 == "T") |
      (A1 == "G" & A2 == "A") |
      (A1 == "G" & A2 == "T")
  )
  message("Remove the ambiguous SNPs..., remaining ", nrow(sumstats), " SNPs.")

  # duplicated SNPs
  sumstats <- dplyr::distinct(sumstats, SNP, .keep_all = TRUE)
  message("Remove the duplicated SNPs..., remaining ", nrow(sumstats), " SNPs.")

  # check estimate of effect size (BETA) and standard error (SE)
  if (beta_col %in% names(sumstats)) {
    names(sumstats)[which(names(sumstats) == beta_col)[1]] <- "BETA"

    if (!is.numeric(sumstats$BETA)) {
      message("Effect size column is not numeric data. Coercing...")
      sumstats$BETA <- as.numeric(as.character(sumstats$BETA))
    }

    sumstats <- dplyr::filter(sumstats, is.finite(BETA))

    # check se
    if (se_col %in% names(sumstats)) {
      names(sumstats)[which(names(sumstats) == se_col)[1]] <- "SE"

      if (!is.numeric(sumstats$SE)) {
        message("Standard error column is not numeric data. Coercing...")
        sumstats$SE <- as.numeric(as.character(sumstats$SE))
      }

      sumstats <- dplyr::filter(sumstats, is.finite(SE) & SE > 0)
    }
  }

  # set BETA as log(OR) when BETA is not available
  if (!beta_col %in% names(sumstats) & or_col %in% names(sumstats)) {
    names(sumstats)[which(names(sumstats) == or_col)[1]] <- "OR"

    message("Infer effect size column from log(OR)...")
    sumstats$BETA <- log(sumstats$OR)

    # check se
    if (se_col %in% names(sumstats)) {
      names(sumstats)[which(names(sumstats) == se_col)[1]] <- "SE"

      if (!is.numeric(sumstats$SE)) {
        message("Standard error column is not numeric data. Coercing...")
        sumstats$SE <- as.numeric(as.character(sumstats$SE))
      }

      sumstats <- dplyr::filter(sumstats, is.finite(SE) & SE > 0)
    }
  }

  # check FREQ
  if (freq_col %in% names(sumstats)) {
    names(sumstats)[which(names(sumstats) == freq_col)[1]] <- "FREQ"
    if (!is.numeric(sumstats$FREQ)) {
      message("Frequency column is not numeric data. Coercing...")
      sumstats$FREQ <- as.numeric(as.character(sumstats$FREQ))
    }

    # remove SNP with allele frequency less than min_freq
    sumstats <- dplyr::filter(sumstats, (FREQ > min_freq & FREQ < (1 - min_freq)) | is.na(FREQ))
    message(
      "Remove the SNPs with allele frequency less than min_freq (keeping NAs)..., remaining ",
      nrow(sumstats), " SNPs."
    )
  }

  # check p-value
  if (p_col %in% names(sumstats) & log_pval) {
    names(sumstats)[which(names(sumstats) == p_col)[1]] <- "P"

    if (!is.numeric(sumstats$P)) {
      message("p-value column is not numeric data. Coercing...")
      sumstats$P <- as.numeric(sumstats$P)
    }

    sumstats$P <- 10^(-sumstats$P)
  } else if (p_col %in% names(sumstats) & (!log_pval)) {
    names(sumstats)[which(names(sumstats) == p_col)[1]] <- "P"

    if (!is.numeric(sumstats$P)) {
      message("p-value column is not numeric data. Coercing...")
      sumstats$P <- as.numeric(sumstats$P)
    }

    sumstats <- dplyr::filter(sumstats, P >= 0 & P <= 1)
    message(
      "Remove the SNPs with p-value < 0 or p-value > 1..., remaining ",
      nrow(sumstats), " SNPs."
    )
  }

  # check z-score
  if (z_col %in% names(sumstats)) {
    names(sumstats)[which(names(sumstats) == z_col)[1]] <- "Z"

    if (!is.numeric(sumstats$Z)) {
      message("z-score column is not numeric data. Coercing...")
      sumstats$Z <- as.numeric(sumstats$Z)
    }

    message("Calculate chi-square column from z-score.")
    sumstats$CHI2 <- (sumstats$Z)^2
  }

  # calculate z-score from p-value
  if ((!"Z" %in% names(sumstats)) & "P" %in% names(sumstats) & "BETA" %in% names(sumstats)) {
    message("Infer z-score from p-value and effect size.")
    sumstats$CHI2 <- stats::qchisq(sumstats$P, 1, lower.tail = FALSE)
    sumstats$Z <- sign(sumstats$BETA) * sqrt(sumstats$CHI2)
  }

  # calculate z-score if no Z, but with BETA and SE
  if ((!"Z" %in% names(sumstats)) &
    ("BETA" %in% names(sumstats) & "SE" %in% names(sumstats))) {
    message("Infer z-score from effect size and standard error.")
    sumstats$Z <- sumstats$BETA / sumstats$SE
    sumstats$CHI2 <- (sumstats$Z)^2
  }

  # calculate chi-square if no CHI2
  if (!"CHI2" %in% names(sumstats)) {
    sumstats$CHI2 <- (sumstats$Z)^2
  }

  # calculate p-value if no P
  if (!"P" %in% names(sumstats)) {
    sumstats$P <- stats::pchisq(sumstats$CHI2, 1, lower.tail = FALSE)
  }

  if (!"Z" %in% names(sumstats)) {
    stop("no information for z-score!")
  }

  # check missing data
  sumstats <- dplyr::select(sumstats, dplyr::one_of(out))
  sumstats <- tidyr::drop_na(sumstats)
  message("Remove missing values..., remaining ", nrow(sumstats), " SNPs.")

  n_min <- mean(sumstats$N) - 5 * stats::sd(sumstats$N)
  n_max <- mean(sumstats$N) + 5 * stats::sd(sumstats$N)
  sumstats <- dplyr::filter(sumstats, N >= n_min & N <= n_max)
  message(
    "Remove SNPs with sample size 5 standard deviations away from the mean..., remaining ",
    nrow(sumstats), " SNPs."
  )

  if (is.null(chi2_max)) {
    chi2_max <- max(c(80, stats::median(sumstats$N) / 1000))
  }
  sumstats <- dplyr::filter(sumstats, CHI2 < chi2_max)
  message(
    "Remove SNPs with chi-square > chi2_max..., remaining ",
    nrow(sumstats), " SNPs."
  )

  message("The formatted data has ", nrow(sumstats), " lines.")

  return(sumstats)
}

#' Get the size of the credible set
#'
#' @param proportions A numeric vector of proportions.
#' @param threshold A numeric value specifying the threshold for the credible set. The default is `0.95`.
#'
#' @return The size of the credible set.
#'
#' @export
getCredibleSetSize <- function(proportions, threshold = 0.95) {
  proportions <- sort(proportions, decreasing = TRUE)
  proportion_sum <- 0
  for (i in seq_along(proportions)) {
    proportion_sum <- proportion_sum + proportions[i]
    if (proportion_sum >= threshold) {
      return(i)
    }
  }
}
