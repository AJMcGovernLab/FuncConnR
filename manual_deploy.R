# Manual Deployment - Alternative Method
# If automated deployment fails, try this step-by-step approach

library(rsconnect)

# Method 1: Set credentials step by step
cat("Manual deployment method:\n")
cat("1. Make sure you have fresh credentials from https://www.shinyapps.io/admin/#/tokens\n")
cat("2. Run each command below one by one:\n\n")

cat("# Step 1: Set account info (replace with your actual values)\n")
cat('rsconnect::setAccountInfo(name="ajmg", token="YOUR_TOKEN", secret="YOUR_SECRET")\n\n')

cat("# Step 2: Check connection\n")
cat("rsconnect::accounts()\n\n")

cat("# Step 3: Deploy app\n")
cat('rsconnect::deployApp(appName="consensusconnectr", forceUpdate=TRUE)\n\n')

# Method 2: Interactive setup
cat("Or try interactive setup:\n")
cat("rsconnect::deployApp()  # This will prompt for credentials\n\n")

# Method 3: Check existing accounts
existing_accounts <- rsconnect::accounts()
if (nrow(existing_accounts) > 0) {
  cat("Existing accounts found:\n")
  print(existing_accounts)
  cat("\nYou can deploy directly with:\n")
  cat('rsconnect::deployApp(appName="consensusconnectr", forceUpdate=TRUE)\n')
} else {
  cat("No existing accounts found. You need to set up authentication first.\n")
}

# Troubleshooting tips
cat("\n=== TROUBLESHOOTING TIPS ===\n")
cat("1. Token and Secret must be DIFFERENT values\n")
cat("2. Make sure you're using the latest credentials from shinyapps.io\n")
cat("3. Try removing old credentials: rsconnect::forgetDeployment()\n")
cat("4. Check your internet connection\n")
cat("5. Try deploying from RStudio GUI instead\n")