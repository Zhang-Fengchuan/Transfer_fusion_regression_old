library(foreign)
library(car)
library(stringr)

bed <- read.csv("C:/Users/张丰川/Desktop/real_data/pancreatic_0_2010-2012_black.csv", sep = ',', header = TRUE)#目标州
#bed <- read.csv("C:/Users/张丰川/Desktop/real_data/pancreatic_1_2010-2012_black.csv", sep = ',', header = TRUE)#加州
#bed <- read.csv("C:/Users/张丰川/Desktop/real_data/pancreatic_2_2010-2012_black.csv", sep = ',', header = TRUE)#乔治亚
#bed <- read.csv("C:/Users/张丰川/Desktop/real_data/pancreatic_3_2010-2012_black.csv", sep = ',', header = TRUE)#路易斯安那
#bed <- read.csv("C:/Users/张丰川/Desktop/real_data/pancreatic_4_2010-2012_black.csv", sep = ',', header = TRUE)#德克萨斯

#View(bed)
head(bed)
Data = cbind(bed$Age.recode.with.single.ages.and.90.,   bed$Sex, bed$Year.of.diagnosis,   bed$Behavior.code.ICD.O.3,
            bed$AYA.site.recode.2020.Revision,   bed$Primary.Site,   bed$Histologic.Type.ICD.O.3,
            bed$Derived.Summary.Grade.2018..2018..,   bed$Grade.Clinical..2018..,   bed$Grade.Pathological..2018..,
            bed$Laterality,   bed$ICD.O.3.Hist.behav,   bed$Combined.Summary.Stage..2004..,
            bed$RX.Summ..Surg.Prim.Site..2003..,  bed$RX.Summ..Surg.Rad.Seq..2006..,   bed$Radiation.recode..2003..,
            bed$Chemotherapy.recode..yes..no.unk...2004..,   bed$RX.Summ..Systemic.Sur.Seq..2007..,
            bed$Time.from.diagnosis.to.treatment.in.days.recode,  bed$Breast.Subtype..2010..,   bed$Derived.HER2.Recode..2010..,
            bed$ER.Status.Recode.Breast.Cancer..2010..,   bed$PR.Status.Recode.Breast.Cancer..2010..,  bed$First.malignant.primary.indicator,
            bed$Total.number.of.in.situ.malignant.tumors.for.patient,  bed$Year.of.death.recode,  bed$Marital.status.at.diagnosis,
            bed$Median.household.income.inflation.adj.to.2022);Data

colnames(Data) <- c("Age.recode.with.single.ages.and.90.", "Sex", "Year.of.diagnosis", "Behavior.code.ICD.O.3",
                    "AYA.site.recode.2020.Revision", "Primary.Site", "Histologic.Type.ICD.O.3",
                    "Derived.Summary.Grade.2018..2018..", "Grade.Clinical..2018..",  "Grade.Pathological..2018..",
                    "Laterality", "ICD.O.3.Hist.behav", "Combined.Summary.Stage..2004..",
                    "Summ..Surg.Prim.Site..2003..", "RX.Summ..Surg.Rad.Seq..2006..", "Radiation.recode..2003..",
                    "Chemotherapy.recode..yes..no.unk...2004..", "RX.Summ..Systemic.Sur.Seq..2007..",
                    "Time.from.diagnosis.to.treatment.in.days.recode", "Breast.Subtype..2010..", "Derived.HER2.Recode..2010..",
                    "Status.Recode.Breast.Cancer..2010..", "PR.Status.Recode.Breast.Cancer..2010..",  "First.malignant.primary.indicator",
                    "Total.number.of.in.situ.malignant.tumors.for.patient", "Year.of.death.recode", "Marital.status.at.diagnosis",
                    "Median.household.income.inflation.adj.to.2022")  # 设置列名

colnames(Data) <- c("age", "sex", "year_of_diagnosis", "behavior", "site_detail", "site_primary", "type", "grade_summary",
                    "grade_clinical", "grade_pathological", "laterality", "hist_behavior", "summary_stage", "surg", "surg_rad_seq", "radiation",
                    "chemotherapy", "systemic_surg_seq", "time_diag_treat", "subtype", "Her2", "ER", "PR", "first_if", "total_number", "year_of_death",
                    "marital_status", "median_income")  # 设置列名
# View(Data)
Data <- data.frame(Data)
names(Data)


# 年龄
Data$age <- gsub(" years", "", Data$age);
Data$age <- as.numeric(Data$age);Data$age
age_mean = mean(Data$age);sd_mean = sd(Data$age);
Data$age <- (Data$age-age_mean)/sd_mean;Data$age


# 性别
Data$sex <- ifelse(Data$sex=="Female",0,ifelse(Data$sex=="Male",1, NA));Data$sex

# 诊断年份
Data$year_of_diagnosis <- as.numeric(Data$year_of_diagnosis);Data$year_of_diagnosis


# 肿瘤性质    Benign良性的   Borderline malignancy潜在恶性的   In situ原位癌；   Malignant恶性的
Data$behavior <- recode(Data$behavior, "'Benign'=0;'Borderline malignancy'=1;'In situ'=2;'Malignant'=3; else = NA");Data$behavior
Data$behavior[which(Data$behavior == 0)] = 0;
Data$behavior[which(Data$behavior != 0)] = 1;

#0 1 1 1 良性为0,非良性为1





#######
# 肿瘤具体位置编码  1:乳腺浸润管 2:乳腺癌 3:乳房 - 小叶 4:乳房 - 叶状体 5:乳房 - 叶状体 
# 肿瘤具体位置编码  6:乳房 - 佩吉特 7:乳腺导管 8:乳房 - 化生 9: 乳房 - 炎症 10:乳房 - 其他
# 肿瘤具体位置编码  数值代表的具体位置查询https://seer.cancer.gov/ayarecode/aya-2020.html   
numbers_only <- regmatches(Data$site_detail, gregexpr("[0-9]+", Data$site_detail));numbers_only
numbers_only <- sapply(numbers_only, paste, collapse = "")
Data$site_detail <- numbers_only;Data$site_detail
Data$site_detail <- as.numeric(Data$site_detail);Data$site_detail
# 如果不是乳腺癌取消这一句
#Data$site_detail <- recode(Data$site_detail, "961=1;962=2;963=3;964=4;965=5;966=6;
#                           967=7;968=8;969=9;9610=10; else = NA");Data$site_detail 
#Data$site_detail[c(which(Data$site_detail == 1),which(Data$site_detail == 2),which(Data$site_detail == 7))] = 1;
#Data$site_detail[setdiff(c(1:length(Data$site_detail)),c(which(Data$site_detail == 1),
#                                                     which(Data$site_detail == 2),which(Data$site_detail == 7)))] = 0;Data$site_detail

# 1 浸润型 癌细胞已扩散到乳腺导管或小叶以外的组织
# 2 非浸润型 癌细胞仍局限在乳腺导管或小叶内





# 胰腺癌具体位置
#1 神经内分泌肿瘤（NET） 2神经内分泌癌（NEC） 3其他神经内分泌肿瘤 4胰腺腺癌 5其他胰腺癌 6其他胃肠道癌症
#0 神经内分泌肿瘤(NETs) 1 2 3           1 非神经内分泌肿瘤(NETs) 4 5 6
Data$site_detail <- recode(Data$site_detail, "11=1;93912=2;93913=3;9392=4;9393=5;9310=6;else = NA");Data$site_detail 
Data$site_detail[setdiff(c(1:length(Data$site_detail)),c(which(Data$site_detail == 4),
                                                         which(Data$site_detail == 5),which(Data$site_detail == 6)))] = 0;Data$site_detail
Data$site_detail[c(which(Data$site_detail == 4),which(Data$site_detail == 5),which(Data$site_detail == 6))] = 1;Data$site_detail












# 肿瘤原发位置 数值代表的具体原发位置查询https://seer.cancer.gov/ayarecode/aya-2020.html 
Data$site_primary <- as.numeric(Data$site_primary);Data$site_primary



# 肿瘤的组织学类型 数值代表的具体组织学类型查询https://seer.cancer.gov/ayarecode/aya-2020.html
Data$type <- as.numeric(Data$type);Data$type


# 肿瘤的总结分级  1:I期 低等级,高分化  2:II期 中等级,中分化  3:III期 高等级 低分化 4: 未分化 9:未知
Data$grade_summary <- as.numeric(Data$grade_summary);Data$grade_summary
Data$grade_summary[which(Data$grade_summary == 1)] = 1;
Data$grade_summary[c(which(Data$grade_summary == 2),which(Data$grade_summary == 3),which(Data$grade_summary == 4))] = 0;
Data$grade_summary[which(Data$grade_summary == 9)] = NA;


# 肿瘤的临床分级  1:I期 低等级,高分化  2:II期 中等级,中分化  3:III期 高等级 低分化 4: 未分化 9:未知
Data$grade_clinical <- as.numeric(Data$grade_clinical);Data$grade_clinical
# 0 1 1 1
Data$grade_clinical[which(Data$grade_clinical == 1)] = 1;
Data$grade_clinical[c(which(Data$grade_clinical == 2),which(Data$grade_clinical == 3),which(Data$grade_clinical == 4))] = 0;
Data$grade_clinical[which(Data$grade_clinical == 9)] = NA;



# 肿瘤的病理分级  1:I期 低等级,高分化  2:II期 中等级,中分化  3:III期 高等级 低分化 4: 未分化 9:未知
Data$grade_pathological <- as.numeric(Data$grade_pathological);Data$grade_pathological
# 0 1 1 1
Data$grade_pathological[which(Data$grade_pathological == 1)] = 1;
Data$grade_pathological[c(which(Data$grade_pathological == 2),which(Data$grade_pathological == 3),which(Data$grade_pathological == 4))] = 0;
Data$grade_pathological[which(Data$grade_pathological == 9)] = NA;

# 肿瘤的位置侧性  0:肿瘤不是成对的位置 1:右侧 2:左侧 3:只有一侧但未指明 4:双侧 5:中线 6:没有信息
Data$laterality <- recode(Data$laterality, "'Not a paired site'=0;'Right - origin of primary'=1;'Left - origin of primary'=2;
                          'Only one side - side unspecified'=3;'Bilateral, single primary'=4; 'Paired side： midline tumor'=5;
                          'Paired site, but no information concerning laterality'=6; else = NA");Data$laterality
#



# 肿瘤的组织学类型 详细参考 See SEER*Stat dictionary for labels
numbers_only <- regmatches(Data$hist_behavior, gregexpr("[0-9]+", Data$hist_behavior));numbers_only
numbers_only <- sapply(numbers_only, paste, collapse = "")
Data$hist_behavior <- numbers_only;Data$hist_behavior
Data$hist_behavior <- as.numeric(Data$hist_behavior);Data$hist_behavior




# 肿瘤扩散程度  0: In situ原位癌 1:Localized局限性 2:Regional区域性 3:Distant远处 9:Unknown/unstaged未知
Data$summary_stage <- recode(Data$summary_stage,"'In situ'=0; 'Localized'=1; 'Regional'=2;
                             'Distant'=3; 'Unknown/unstaged'=9; else = NA");Data$summary_stage
# 0 0 0 1
Data$summary_stage[which(Data$summary_stage != 3)] = 0;
Data$summary_stage[which(Data$summary_stage == 3)] = 1;





# 原发部位手术切除情况 0(0/00):没有手术 1(10-19):局部肿瘤破坏 2(20-29):局部肿瘤切除 3(30-39):部分或次全切除 
# 原发部位手术切除情况 4(40-49):全切或近全切除 5(50-59):包括连续切除的手术  6(60-69):保肢手术 7(70-79)截肢手术
#原发部位手术切除情况  8(80-89):根治手术  9(90):手术未另行说明 10(99):未知是否进行手术
Data$surg <- as.numeric(Data$surg);Data$surg
Data$surg <- ifelse(Data$surg==0,0,ifelse(10<=Data$surg&Data$surg<=19,1,ifelse(20<=Data$surg&Data$surg<=29,2,
ifelse(30<=Data$surg&Data$surg<=39,3,ifelse(40<=Data$surg&Data$surg<=49,4,ifelse(50<=Data$surg&Data$surg<=59,5,
ifelse(60<=Data$surg&Data$surg<=69,6,ifelse(70<=Data$surg&Data$surg<=79,7,ifelse(80<=Data$surg&Data$surg<=89,8,
ifelse(Data$surg&Data$surg==90,9,ifelse(Data$surg&Data$surg==99,10,NA)))))))))));Data$surg
Data$surg[which(Data$surg < 1)] = 0;
Data$surg[which(Data$surg >= 1)] = 1;
#0 0 0 0 1 1 1 1 1 1 1



# 手术和放疗的顺序信息 0:'No radiation and/or no surgery; unknown if surgery and/or radiation given'没有进行放疗和/或癌症导向手术
# 1:'Radiation before surgery' 手术前放射治疗  2:'Radiation after surgery' 手术后放射治疗
# 3:'Radiation before and after surgery' 手术前后放射治疗  4:'Intraoperative radiation' 术中放射治疗
# 5: Intraoperative radiation with other radiation given before/after surgery 术中放射治疗并在手术前后进行其他放射治疗
# 6: Surgery both before and after radiation 放射治疗前后的手术
# 7: Unknown if surgery and/or radiation given (未知是否进行了手术和/或放射治疗)
sting_0 <-  'No radiation and/or no surgery; unknown if surgery and/or radiation given';sting_0
Data$surg_rad_seq <- recode(Data$surg_rad_seq, "sting_0 = 0;
                            'Radiation prior to surgery'=1; 'Radiation after surgery'=2; 'Radiation before and after surgery'=3;
                            'Intraoperative radiation'=4; 'Intraoperative radiation with other radiation given before/after surgery'=5;
                            'Surgery both before and after radiation'=6; 'Sequence unknown, but both were given'=7; else = NA");Data$surg_rad_seq
Data$surg_rad_seq[which(Data$surg_rad_seq == 0)] = 0;
Data$surg_rad_seq[which(Data$surg_rad_seq != 0)] = 1;
#0  else 1
#0 没有同时接受手术和放疗   1 同时接受手术和放疗



# 放射治疗信息 'None/Unknown':0 未进行放疗 'Beam radiation':1 光束辐射
# 'Radioactive implants (includes brachytherapy) (1988+)':2 放射性植入物 'Radioisotopes (1988+)':3 放射性同位素
# 'Combination of beam with implants or isotopes':4 光束放疗与植入物或同位素的组合
# 'Radiation, NOS  method or source not specified':5 未特指的放疗 
# 'Refused (1988+)':6 拒绝放疗   'Recommended, unknown if administered':7 未知是否进行了放疗
Data$radiation <- recode(Data$radiation, "'None/Unknown'=0;'Beam radiation'=1;'Radioactive implants (includes brachytherapy) (1988+)'=2;
                         'Radioisotopes (1988+)'=3;'Combination of beam with implants or isotopes'=4;'Radiation, NOS  method or source not specified'=5;
                         'Refused (1988+)'=6; 'Recommended, unknown if administered'=7; else = NA");Data$radiation

#0 0 0 0 1 1 1 1 
#0 未接受放疗或未知  1 接受了某种放疗
Data$radiation[c(which(Data$radiation == 0),which(Data$radiation == 6),which(Data$radiation == 7))] = 0;
Data$radiation[setdiff(c(1:length(Data$radiation)),c(which(Data$radiation == 0),
                                                     which(Data$radiation == 6),which(Data$radiation == 7)))] = 1;Data$radiation





# 是否接受化疗 'No/Unknown':0 未接受化疗  'Yes':1 接受化疗
Data$chemotherapy <- recode(Data$chemotherapy, "'No/Unknown'=0; 'Yes'=1; else = NA");Data$chemotherapy





# 系统性治疗和手术的顺序  'No systemic therapy and/or surgical procedures':0 无系统性治疗和/或手术
# 'Systemic therapy before surgery':1 系统性治疗在手术前 
# 'Systemic therapy after surgery':2 系统性治疗在手术后
# 'Systemic therapy both before and after surgery':3 系统性治疗在手术前后 
# 'Intraoperative systemic therapy'4 术中系统性治疗
# 'Intraop systemic rx & oth systemic rx before/after surg':5 术中系统性治疗和/或手术前后其他系统性治疗
# 'Surgery both before and after systemic therapy':6 系统性治疗前后手术
# 'Sequence unknown':7 顺序未知 
# 'Blank(s)':8 空白
Data$systemic_surg_seq <- recode(Data$systemic_surg_seq, "'No systemic therapy and/or surgical procedures'=0;
                                 'Systemic therapy before surgery'=1; 'Systemic therapy after surgery'=2; 
                                 'Systemic therapy both before and after surgery'=3; 'Intraoperative systemic therapy'=4;
                                 'Intraop systemic rx & oth systemic rx before/after surg'=5;
                                 'Surgery both before and after systemic therapy'=6;'Sequence unknown'=7;'Blank(s)'=8; else = NA");Data$systemic_surg_seq
# 0 else 1
# 0 没有进行系统性治疗 1 进行了系统性治疗
Data$systemic_surg_seq[c(which(Data$systemic_surg_seq == 0),which(Data$systemic_surg_seq == 7),which(Data$systemic_surg_seq == 8))] = 0;
Data$systemic_surg_seq[setdiff(c(1:length(Data$systemic_surg_seq)),c(which(Data$systemic_surg_seq == 0),
                              which(Data$systemic_surg_seq == 7),which(Data$systemic_surg_seq == 8)))] = 1;Data$systemic_surg_seq





# 从诊断到治疗的时间
Data$time_diag_treat <- as.numeric(Data$time_diag_treat);Data$time_diag_treat
#mean_value <- mean(Data$time_diag_treat, na.rm = TRUE)
# 用均值填补NA
#Data$time_diag_treat <- ifelse(is.na(Data$time_diag_treat), mean_value, Data$time_diag_treat);Data$time_diag_treat




# 乳腺癌亚型 'HR+/HER2+':1   'HR+/HER2-':2  'HR-/HER2+':3  HR-/HER2-':4  'Unknown':5
Data$subtype <- recode(Data$subtype, "'HR+/HER2+'=1; 'HR+/HER2-'= 2; 'HR-/HER2+'=3; 'HR-/HER2-'=4; 'Unknown'=5; else = NA");Data$subtype




# HER2 人类表皮生长因子受体2  'Positive':1 阳性  'Negative':2 阴性 'Borderline/Unknown': 边界或未知
Data$Her2 <- recode(Data$Her2, "'Positive'=1; 'Negative'=2; 'Borderline/Unknown'=3; else = NA");Data$Her2



# ER 雌激素受体 'Positive':1 阳性  'Negative':2 阴性 'Borderline/Unknown': 边界或未知
Data$ER <- recode(Data$ER, "'Positive'=1; 'Negative'=2; 'Borderline/Unknown'=3; else = NA");Data$ER



# PR 孕激素受体 'Positive':1 阳性  'Negative':2 阴性 'Borderline/Unknown': 边界或未知
Data$PR <- recode(Data$PR, "'Positive'=1; 'Negative'=2; 'Borderline/Unknown'=3; else = NA");Data$PR




# 是否为首个原发癌 'No':0  'Yes':1
Data$first_if <- recode(Data$first_if, "'No'=0; 'Yes'=1; else = NA");Data$first_if



# 每个病例诊断出的原位和恶性肿瘤的总次数
Data$total_number <- as.numeric(Data$total_number);Data$total_number





# 诊断时婚姻状态 'Single (never married)':0 单身     'Married (including common law)':1 已婚
# 诊断时婚姻状态 'Separated':2 分居                 'Divorced':3 离婚
# 诊断时婚姻状态 'Widowed':4 丧偶                  'Unmarried or Domestic Partner':5 未婚或同居
# 诊断时婚姻状态 'Unknown' 6 未知
Data$marital_status <- recode(Data$marital_status,"'Single (never married)' = 0; 'Married (including common law)' = 1;
                              'Separated'=2; 'Divorced'=3; 'Widowed'=4; 'Unmarried or Domestic Partner'=5;'Unknown'=6; else = NA");Data$marital_status
Data$marital_status[which(Data$marital_status == 1)] = 1;
Data$marital_status[which(Data$marital_status != 1)] = 0;
# one or two






# 家庭收入中位数
for (i in 1:length(Data$median_income)) {
  string_no_commas <- gsub(",", "", Data$median_income[i]); string_no_commas
  numbers <- as.numeric(unlist(regmatches(string_no_commas, gregexpr("[0-9]+", string_no_commas))));numbers
  Data$median_income[i] <- as.numeric(mean(numbers));Data$median_income[i]
}
Data$median_income <- as.numeric(Data$median_income);Data$median_income
mean <- mean(Data$median_income);
sd <- sd(Data$median_income);
Data$median_income <- (Data$median_income-mean)/sd;Data$median_income







# 每一列进行缺失值中位数填补
# 定义一个函数，用于将NA替换为列的中位数
fill_na_with_median <- function(x) {
  median_value <- median(x, na.rm = TRUE)  # 计算中位数，忽略NA
  x[is.na(x)] <- median_value  # 将NA替换为中位数
  return(x)
}

# 使用lapply对数据框的每一列应用该函数
Data <- as.data.frame(lapply(Data, fill_na_with_median))
# 对发现到治疗时间进行标准化
mean <- mean(Data$time_diag_treat);
sd <- sd(Data$time_diag_treat);
Data$time_diag_treat <- (Data$time_diag_treat-mean)/sd;Data$time_diag_treat
mean <- mean(Data$total_number);
sd <- sd(Data$total_number);
Data$total_number <- (Data$total_number-mean)/sd;Data$total_number





# 计算生存时间
# 死亡年份 'Alive at last contact':2021+  
Data$year_of_death <- recode(Data$year_of_death, "'Alive at last contact'=9999");Data$year_of_death
Data$year_of_death <- as.numeric(Data$year_of_death);Data$year_of_death
Data$year_of_diagnosis
Data$year_of_death
time <- matrix(0, nrow = length(Data$year_of_death), ncol = 1);time
delta <- matrix(0, nrow = length(Data$year_of_death), ncol = 1);delta
for (i in c(1:length(Data$year_of_death))) {
  if (Data$year_of_death[i]== 9999){
    delta[i,1] <- 0;
    time[i,1] <- log(2021 - Data$year_of_diagnosis[i]);
  }else{
    delta[i,1] <- 1;
    if ((Data$year_of_death[i] == Data$year_of_diagnosis[i])){
      time[i,1] <- log(0.5);
    }else{
      time[i,1] <- log(Data$year_of_death[i] - Data$year_of_diagnosis[i]);
    }
  }
}
#time #delta #Data
n = length(time);n#样本量
sorted <- order(time);sorted
Data <- Data[sorted, ];#View(Data)
delta <- delta[sorted,];delta
weight <- matrix(0, nrow = n, ncol = 1);weight
for (i in c(1:n)) {
  if(i == 1){
    weight[i,1] <- delta[i]/n;
  }else{
    weight_change <- 1;
    for (j in c(1:(i-1))) {
      weight_change <- weight_change*((n-j)/(n-j+1))^(delta[j]);
    }
    weight_change <- weight_change*delta[i]/(n-i+1);
    weight[i,1] <- weight_change;
  }
}
# weight
time_mean <- sum(weight*time)/sum(weight);time_mean
time <- sqrt(weight)*(time-time_mean);time
Y <- time;#Y
# Data$laterality
X <- cbind.data.frame(Data$age, Data$sex, Data$behavior, Data$site_detail,
                      Data$summary_stage, Data$surg, Data$surg_rad_seq,  Data$radiation,
                      Data$chemotherapy, Data$systemic_surg_seq, Data$time_diag_treat,
                      Data$first_if, Data$total_number,  Data$marital_status, 
                      Data$median_income);#View(X)
X2 <- cbind.data.frame(Data$site_primary, Data$type, Data$grade_summary, Data$grade_clinical,
                       Data$grade_pathological, Data$hist_behavior, Data$subtype, Data$Her2,
                       Data$ER, Data$PR, Data$grade_summary);#View(X2)
X_mean <- matrix(0, nrow = 1, ncol = ncol(X));X_mean
for (i in c(1:nrow(X))) {
  X_mean <- X_mean + X[i, ]*weight[i, 1];
}
X_mean <- X_mean/sum(weight);X_mean
for (i in c(1:nrow(X))) {
  for (j in c(1:ncol(X))) {
    X[i, j] <- sqrt(weight[i, ])*(X[i, j] - X_mean[ ,j]);
  }
};#View(X)




Data_out <- data.frame(
  Y, X, X2
);#View(Data_out)
write.csv(Data_out, "C:/Users/张丰川/Desktop/real_data/Data_target.csv", row.names = FALSE)


#Data_out <- data.frame(
#  Y, X, X2
#);#View(Data_out)
#write.csv(Data_out, "C:/Users/张丰川/Desktop/real_data/Data_1.csv", row.names = FALSE)


#Data_out <- data.frame(
#  Y, X, X2
#);#View(Data_out)
#write.csv(Data_out, "C:/Users/张丰川/Desktop/real_data/Data_2.csv", row.names = FALSE)


#Data_out <- data.frame(
#  Y, X, X2
#);#View(Data_out)
#write.csv(Data_out, "C:/Users/张丰川/Desktop/real_data/Data_3.csv", row.names = FALSE)


#Data_out <- data.frame(
#  Y, X, X2
#);#View(Data_out)
#write.csv(Data_out, "C:/Users/张丰川/Desktop/real_data/Data_4.csv", row.names = FALSE)
#View(X)












