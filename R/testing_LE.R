data <- read.csv("Data/Data_Félicien/BCI_v5.1.csv", stringsAsFactors = FALSE)

data$date <- as.POSIXct(data$date, format = "%Y-%m-%d %H:%M")
data$hour <- as.numeric(format(data$date, "%H"))
night_data <- subset(data, hour >= 18 | hour < 6)

library(ggplot2)

p = ggplot(night_data, aes(x = date, y = Rnet)) +
  geom_line() +
  theme_minimal() +
  labs(title = "LE over time",
       x = "Date",
       y = "LE")

print(p)
