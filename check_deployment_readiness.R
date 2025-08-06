# ConsensusConnectR - Deployment Readiness Check
# Run this script to verify your app is ready for deployment

cat("=============================================================================\n")
cat("ConsensusConnectR - Deployment Readiness Check\n")
cat("=============================================================================\n\n")

# Check 1: Required packages
cat("1. Checking required packages...\n")
required_packages <- c(
  "shiny", "shinydashboard", "DT", "shinyjs", "colourpicker",
  "ggplot2", "scales", "mice", "igraph", "corrplot", "viridis",
  "psych", "corpcor", "openxlsx", "jsonlite", "rsconnect"
)

missing_packages <- c()
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
    cat("   ❌", pkg, "- NOT INSTALLED\n")
  } else {
    cat("   ✅", pkg, "- OK\n")
  }
}

if (length(missing_packages) > 0) {
  cat("\n   INSTALL MISSING PACKAGES:\n")
  cat("   install.packages(c(", paste0('"', missing_packages, '"', collapse = ", "), "))\n\n")
} else {
  cat("   All packages are installed! ✅\n\n")
}

# Check 2: Essential files
cat("2. Checking essential files...\n")
essential_files <- c("app.R", "modules/analysis_functions.R", "modules/visualization_functions.R", 
                    "modules/ui_components.R", "modules/comprehensive_download_handler.R")

all_files_present <- TRUE
for (file in essential_files) {
  if (file.exists(file)) {
    cat("   ✅", file, "- OK\n")
  } else {
    cat("   ❌", file, "- MISSING\n")
    all_files_present <- FALSE
  }
}

if (!all_files_present) {
  cat("   Some essential files are missing!\n\n")
} else {
  cat("   All essential files present! ✅\n\n")
}

# Check 3: App can load
cat("3. Testing app loading...\n")
tryCatch({
  # Try to source the main app file
  temp_env <- new.env()
  source("app.R", local = temp_env)
  cat("   ✅ App loads successfully!\n\n")
}, error = function(e) {
  cat("   ❌ App failed to load:", e$message, "\n\n")
})

# Check 4: File sizes
cat("4. Checking app size...\n")
app_size <- sum(file.size(list.files(recursive = TRUE, full.names = TRUE), na.rm = TRUE))
app_size_mb <- round(app_size / 1024 / 1024, 2)
cat("   App size:", app_size_mb, "MB\n")

if (app_size_mb > 100) {
  cat("   ⚠️  Large app size - consider excluding unnecessary files\n")
  cat("   Create .shinyappsignore file to exclude large/unnecessary files\n\n")
} else {
  cat("   ✅ App size is reasonable for deployment\n\n")
}

# Summary
cat("=============================================================================\n")
cat("DEPLOYMENT READINESS SUMMARY\n")
cat("=============================================================================\n")

readiness_score <- 0
if (length(missing_packages) == 0) readiness_score <- readiness_score + 1
if (all_files_present) readiness_score <- readiness_score + 1
if (app_size_mb <= 100) readiness_score <- readiness_score + 1

if (readiness_score >= 2) {
  cat("🚀 Your app is ready for deployment!\n\n")
  cat("Next steps:\n")
  cat("1. Get your shinyapps.io credentials from: https://www.shinyapps.io/admin/#/tokens\n")
  cat("2. Edit deploy_app.R with your credentials\n")
  cat("3. Run deploy_app.R to deploy your app\n")
} else {
  cat("⚠️  Please address the issues above before deploying.\n")
}

cat("\n=============================================================================\n")