library('diptest')

# Load and view data
word_report <- read.csv("V:/users/ma09/Projects/AVSpeechMEG/data/derivatives/behavioural/group/word_report.csv")
View(word_report)

# Extract relevant conditions
acc_AOhigh <- word_report$accuracy[(word_report$stim_modality %in% 'AO') & (word_report$acc_clarity %in% 'high')]
acc_AVlow <- word_report$accuracy[(word_report$stim_modality %in% 'AV') & (word_report$acc_clarity %in% 'low')]
acc_VO <- word_report$accuracy[(word_report$stim_modality %in% 'VO')]

# Plot (AV_low-AO_high) vs. VO for sanity check
plot(acc_VO, acc_AVlow - acc_AOhigh)

# Perform dip test on (AV_low-AO_high)
dip.test(acc_AVlow-acc_AOhigh)

