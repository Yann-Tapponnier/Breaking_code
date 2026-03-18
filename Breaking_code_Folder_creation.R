# Path to the folder you want to check/create
folder <- "path/to/your_folder"

# Check if the folder exists; if not, create it
if (!dir.exists(folder)) {
  dir.create(folder, recursive = TRUE)
  message("Folder created: ", folder)
} else {
  message("Folder already exists: ", folder)
}
