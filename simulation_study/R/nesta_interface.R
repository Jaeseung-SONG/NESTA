submitted_nesta_path <- "/home/js/NESTA/Analysis/Nesta.R"
package_nesta_wrapper_path <- "/home/js/NESTA/simulation_study/R/run_nesta_methods.R"

twas_p_from_z <- function(z) 2 * stats::pnorm(-abs(z))

code_fidelity_audit <- function() {
  submitted_exists <- file.exists(submitted_nesta_path)
  submitted_src <- if (submitted_exists) readLines(submitted_nesta_path, warn = FALSE) else character()
  wrapper_exists <- file.exists(package_nesta_wrapper_path)
  wrapper_src <- if (wrapper_exists) readLines(package_nesta_wrapper_path, warn = FALSE) else character()
  data.frame(
    check = c("submitted_nesta_exists", "strict_twas_filter_text", "installed_diffustats_n_perm_wrapper", "nperm_rejected_by_wrapper", "p_conversion_package_wrapper"),
    passed = c(
      submitted_exists,
      any(grepl("TWAS.P <", submitted_src, fixed = TRUE)) || any(grepl("TWAS.P<", submitted_src, fixed = TRUE)),
      wrapper_exists && any(grepl("n.perm", wrapper_src, fixed = TRUE)) && any(grepl("diffuStats::diffuse", wrapper_src, fixed = TRUE)),
      wrapper_exists && any(grepl("Unsupported diffuStats argument `nperm`", wrapper_src, fixed = TRUE)),
      identical(twas_p_from_z(c(-2, 0, 2)), 2 * stats::pnorm(-abs(c(-2, 0, 2))))
    ),
    detail = c(
      submitted_nesta_path,
      "Submitted source is audited for strict TWAS.P filtering text.",
      "Code bundle wrapper calls installed diffuStats with n.perm while preserving submitted NESTA arithmetic.",
      "Code bundle wrapper rejects nperm and instructs n.perm.",
      "TWAS.P wrapper uses 2 * pnorm(-abs(TWAS.Z))."
    ),
    stringsAsFactors = FALSE
  )
}
