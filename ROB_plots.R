library(robvis)

# Read in data
robins_data = read.csv("robins_study_data.csv", header = TRUE, sep = ';')
# This is just the same data with abbreviated names to display in traffic plot
robins_data2 = read.csv("robins_study_data2.csv", header = TRUE, sep = ';')

# make weighting of studies equal for all studies
robins_data$Weight = rep(1, 18)
robins_data2$Weight = rep(1, 18)

# needed to shorted this one abbreviation to fit on plot
robins_data2[7,1] <- "HGN"

# Summary Plot
rob_summary(data = robins_data, 
            tool = "ROBINS-I")

# Traffic Light Plot
rob_traffic_light(data = robins_data2, 
                  tool = "ROBINS-I", psize = 8)

