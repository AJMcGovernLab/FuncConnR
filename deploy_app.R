# ConsensusConnectR Deployment Script for shinyapps.io
# Run this script to deploy your app to shinyapps.io

# Install rsconnect if not already installed
if (!requireNamespace("rsconnect", quietly = TRUE)) {
  install.packages("rsconnect")
}

library(rsconnect)

# =============================================================================
# STEP 1: CONFIGURE YOUR SHINYAPPS.IO CREDENTIALS
# =============================================================================
# Get these from: https://www.shinyapps.io/admin/#/tokens

# REPLACE THESE WITH YOUR ACTUAL CREDENTIALS:
SHINYAPPS_USERNAME <- "ajmg"      # Your shinyapps.io username
SHINYAPPS_TOKEN <- "EA4421F38D8798A0008032DC23C956B2"            # Token from shinyapps.io
SHINYAPPS_SECRET <- "XBdu/ctJ8t8d992bupoZF22s3uT6IqwwxcGEYzzE"          # Secret from shinyapps.io

# Set up authentication (run this once)
cat("Setting up shinyapps.io authentication...\n")
rsconnect::setAccountInfo(
  name = SHINYAPPS_USERNAME,
  token = SHINYAPPS_TOKEN,
  secret = SHINYAPPS_SECRET
)

# =============================================================================
# STEP 2: CHECK REQUIRED PACKAGES
# =============================================================================
cat("Checking required packages...\n")

required_packages <- c(
  "shiny", "shinydashboard", "DT", "shinyjs", "colourpicker",
  "ggplot2", "scales", "mice", "igraph", "corrplot", "viridis",
  "psych", "corpcor", "openxlsx", "jsonlite"
)

missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages)
} else {
  cat("All required packages are installed.\n")
}

# =============================================================================
# STEP 3: PREPARE APP FOR DEPLOYMENT
# =============================================================================
cat("Preparing app for deployment...\n")

# Create .shinyappsignore file to exclude unnecessary files
ignore_content <- c(
  "*.log",
  "*.Rhistory",
  ".RData", 
  ".Ruserdata",
  "packrat/",
  "rsconnect/",
  "DEPLOYMENT_GUIDE.md",
  "deploy_app.R",
  "*.md",
  "Rplots.pdf"
)

writeLines(ignore_content, ".shinyappsignore")
cat("Created .shinyappsignore file\n")

# =============================================================================
# STEP 4: DEPLOY THE APP
# =============================================================================
cat("Deploying app to shinyapps.io...\n")
cat("This may take 5-15 minutes depending on your internet connection.\n")

# Deploy the app
rsconnect::deployApp(
  appName = "consensusconnectr",                                    # App name on shinyapps.io
  appTitle = "ConsensusConnectR - Multimethod Consensus Analysis", # Display title
  account = SHINYAPPS_USERNAME,                                     # Your username
  forceUpdate = TRUE,                                               # Force update if exists
  launch.browser = TRUE                                             # Open in browser when done
)

cat("\n=============================================================================\n")
cat("DEPLOYMENT COMPLETE!\n")
cat("=============================================================================\n")
cat("Your app should now be available at:\n")
cat(paste0("https://", SHINYAPPS_USERNAME, ".shinyapps.io/consensusconnectr/\n"))
cat("\nIf deployment failed, check the error messages above and:\n")
cat("1. Verify your credentials are correct\n")
cat("2. Ensure all required packages are installed\n")
cat("3. Check your internet connection\n")
cat("4. Consider upgrading your shinyapps.io plan for larger apps\n")