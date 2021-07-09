library(ggplot2)

# EMT violin plot (SNUH)
df <- data.frame(TYPE=c("WT", "BRCA", "BRCA", "WT", "BRCA", "WT", "WT", "BRCA", "BRCA", "BRCA", "BRCA", "BRCA", "BRCA", "BRCA", "BRCA", "BRCA", "BRCA", "WT", "BRCA", "BRCA"), EMT_index = c(14.1123922142792, 6.76436432443856, 7.57527263485835, 15.0423739828603, 6.23908032722261, 15.104632507528, 19.270741407756, 10.8424486624009, 10.1941890739711, 6.18055110054678, 9.57532921782817, 5.23457009697109, 6.24503649280527, 7.24717414722779, 6.81105461869981, 10.6080775167307, 9.53260081123286, 12.0776818732693, 9.9763257869376, 4.96275557036637))

ggplot(df, aes(x=TYPE, y=EMT_index, fill=TYPE)) + geom_violin(trim=FALSE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000'))  + theme_classic()

# EMT-related violin plot (SUN-HGOC)
dat2 <- read.csv("CDH1_VIM_TGFB1.csv", row.names=1)
# CDH1
df <- data.frame(Type=c("E", "E", "E", "E", "E", "E", "E", "EMT", "E", "E", "E", "E", "E", "E", "E", "E", "EMT", "EMT", "EMT", "EMT"), CDH1_expression=c(13.9411427147277,	15.506352381592,	14.3413445720896,	13.2327938639303,	13.5166952124584,	15.2891732811763,	13.8803932448834,	11.3158646161657,	15.2046765463172,	14.3627592375224,	15.5611932236906,	15.2467867498199,	14.1546059702608,	14.7990189346941,	14.0236616704698,	12.5891156695508,	12.768656850484,	12.497597383967,	10.9187855696642,	11.0135315054797))
ggplot(df, aes(x=Type, y=CDH1_expression, fill=Type)) + geom_violin(trim=FALSE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000')) + theme_classic()

# VIM
df <- data.frame(Type=c("E", "E", "E", "E", "E", "E", "E", "EMT", "E", "E", "E", "E", "E", "E", "E", "E", "EMT", "EMT", "EMT", "EMT"), VIM_expression=c(15.7716815812217,	15.0167142772187,	14.5038706122843,	16.3140694485307,	15.6691790805546,	15.689635021249,	15.2819413071924,	18.5912624143979,	15.5896060791005,	13.8450965535524,	15.1818615157471,	14.4139545891511,	13.6436198602239,	16.5834032952149,	13.8389635102752,	14.6289069319389,	17.3833431180501,	17.5776614778827,	18.3826067683646,	18.563750891245))
ggplot(df, aes(x=Type, y=VIM_expression, fill=Type)) + geom_violin(trim=FALSE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000')) + theme_classic()

# TGFB1
df <- data.frame(Type=c("E", "E", "E", "E", "E", "E", "E", "EMT", "E", "E", "E", "E", "E", "E", "E", "E", "EMT", "EMT", "EMT", "EMT"), TGFB1_expression=c(9.33839926598193,	9.50610967584407,	9.66166096495095,	9.45121599118046,	10.2359274403974,	9.78774638254001,	9.34218800016942,	10.3631706642984,	9.89270880665189,	9.25944454990826,	9.49125554122715,	8.70428056377323,	9.15424479341798,	10.3433816082702,	8.95421003487029,	10.0737904092824,	10.8462141613539,	11.2577479262528,	10.4713149595074,	11.019761724564))
ggplot(df, aes(x=Type, y=TGFB1_expression, fill=Type)) + geom_violin(trim=FALSE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000')) + theme_classic()
