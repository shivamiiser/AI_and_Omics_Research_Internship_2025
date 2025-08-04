#Working Directory

getwd()


# Inside the project directory, create the following subfolders using R code:
# raw_data, clean_data, scripts, results or Tasks, plots etc

dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("Tasks")
dir.create("plots")


# load the dataset into your R environment
data <- read.csv("D:\\AI_Omics_Internship_2025\\Dataset\\patient_info.csv")

View(data)

# Inspect the structure of the dataset using appropriate R functions
str(data)
summary(data)


# Identify variables with incorrect data types

#gender,  diagnosis,  these attributes can be factor 
# and smoker can be in binary factor as a better data type
# diagnosis can also be represented in binary factor form 1 for Cancer 
# and 0 for normal



# Convert variables to appropriate data types where needed

data$gender <- as.factor(data$gender)
data$diagnosis <- as.factor(data$diagnosis)

str(data)



# Create a new variable for smoking status as a binary factor:

# 1 for "Yes", 0 for "No"


data$smoker_fac <- as.factor(data$smoker)
#numeric number
data$smoker_num <- ifelse(data$smoker_fac == "Yes", 1, 0)
View(data)
data$smoker_num <- as.factor(data$smoker_num) 

df_new = data[,-7]
data = data[,-7]
data = data[,-7]
str(data)






