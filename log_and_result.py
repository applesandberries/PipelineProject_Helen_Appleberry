import os

# Replace 'YourFirstName' and 'YourLastName' with your actual first and last name
directory_name = f"PipelineProject_Helen_Appleberry"
os.system(f"mkdir {directory_name}")
os.chdir(directory_name)

# Your existing automation code here...

# Example: Write log information to the log file in the current directory
with open("PipelineProject.log", "w") as log_file:
    log_file.write("Log")
