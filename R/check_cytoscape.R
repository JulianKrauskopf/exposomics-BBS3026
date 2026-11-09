# Exposomics course — Cytoscape setup check
#
# Run this ONCE in the RStudio Console after installing Cytoscape:
#
# source("R/check_cytoscape.R")
#
# Run this check during your initial course software setup.
# Keep Cytoscape Desktop open while this check runs.
# If either required Cytoscape app is missing, this script installs it.

cat("\n=== Exposomics: Cytoscape setup check ===\n")

if (!requireNamespace("RCy3", quietly = TRUE)) {
  stop(
    "RCy3 is not installed. First run ",
    "source(\"R/setup_packages.R\") in the RStudio Console, ",
    "then run this check again.",
    call. = FALSE
  )
}

cat("\nStep 1 of 3: Connecting RStudio to Cytoscape...\n")

connection <- tryCatch(
  RCy3::cytoscapePing(),
  error = function(e) e
)

if (inherits(connection, "error") ||
    identical(connection, FALSE)) {
  stop(
    "RStudio cannot connect to Cytoscape.\n",
    "1. Start the Cytoscape Desktop application.\n",
    "2. Wait until its main window is fully open.\n",
    "3. Run source(\"R/check_cytoscape.R\") again.\n",
    "If this still fails, check local firewall/security prompts.\n",
    if (inherits(connection, "error")) {
      paste0("Details: ", conditionMessage(connection))
    } else {
      ""
    },
    call. = FALSE
  )
}

cat("OK — RStudio is connected to Cytoscape.\n")

is_installed <- function(app_name) {
  installed <- RCy3::getInstalledApps()
  any(grepl(
    app_name,
    as.character(unlist(installed, use.names = FALSE)),
    ignore.case = TRUE,
    fixed = TRUE
  ))
}

required_apps <- c("WikiPathways", "stringApp")

for (i in seq_along(required_apps)) {
  app <- required_apps[i]
  cat(
    "\nStep ", i + 1L, " of 3: Checking ", app, "...\n",
    sep = ""
  )

  available <- tryCatch(
    is_installed(app),
    error = function(e) {
      stop(
        "Cannot check Cytoscape apps: ",
        conditionMessage(e),
        "\nPlease update Cytoscape or RCy3 and retry.",
        call. = FALSE
      )
    }
  )

  if (available) {
    cat("OK — ", app, " is installed.\n", sep = "")
  } else {
    cat(
      app, " is missing. Installing it now ",
      "(internet required)...\n",
      sep = ""
    )

    tryCatch(
      RCy3::installApp(app),
      error = function(e) {
        stop(
          "Could not install ", app, ": ",
          conditionMessage(e),
          "\nCheck your internet connection, then retry.",
          call. = FALSE
        )
      }
    )

    if (!is_installed(app)) {
      stop(
        app, " was not detected after installation.\n",
        "Restart Cytoscape if necessary, then run the ",
        "setup check again.",
        call. = FALSE
      )
    }

    cat("OK — ", app, " installed successfully.\n", sep = "")
  }
}

cat(
  "\nSUCCESS: Cytoscape and both course apps are ready!\n",
  "You can close Cytoscape now. Reopen it when you ",
  "reach Session 6.\n",
  "\n",
  sep = ""
)
