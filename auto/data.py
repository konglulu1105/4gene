# -*- coding: UTF-8 -*-
import os
import sys
import datetime
import smtplib
from email.mime.text import MIMEText
from email.header import Header

###########################################################################################
UPLOAD1 = "/data01/project/sequence/share_local/mgiseq/10030300210002/Info/Upload"
UPLOAD2 = "/data01/project/sequence/share_local/mgiseq/13030400210051/Info/Upload"
UPLOAD3 = "/data01/project/sequence/share_local/mgiseq/1303400210053/Info/Upload"
UPLOAD4 = "/data01/project/sequence/share_local/mgiseq/130340220038/Info/Upload"
UPLOAD5 = "/data01/project/sequence/share_local/mgiseq/130340220036/Info/Upload"
Sequences = "/data01/project/sequence/2023_data/"
#Sequences = "/data01/home/robot02/Project/sequence/"

TMP = "/data01/home/robot02/.TMP/"
data = "/data01/home/robot02/.TMP/data.log"
sample_dir = "/data01/home/robot02/SampleList/"

UnAnalysis = "/data01/home/robot02/.TMP/UnAnalysis.log"
TotalTrans = "/data01/home/robot02/.TMP/TotalTrans.log"
############################################################################################
def diff(UnAnalysis,TotalTrans):
    os.system("ls %s > %s"%(sample_dir,TotalTrans))
    file_1 = open(TotalTrans,"r")
    file_2 = open(UnAnalysis,"r")
    a = file_1.read().splitlines()
    b = file_2.read().splitlines()
    file_1.close()
    file_2.close()
    ready = ([x for x in a if x not in b])
    if len(ready) == 0:
        raise SystemExit('没有待传输!')
    else:
        return ready
###############################读送样信息单的项目号和芯片号######################################
def load_sample(sample_list):
    with open(sample_list, "r") as list:
        for line in list:
            line = line.split("\t")
            project = line[0]
            if len(line) >= 15:
                chip = line[14].strip("\n")
            else:
                chip = "Not"
    return project, chip


###############################列出三台机器下机数据的log文件#####################################
def Log_compare(UPLOAD, chip):
    os.system("ls %s > %s" % (UPLOAD, data))
    with open(data, "r") as file:
        for line in file:
            if chip in line:
                PATH = UPLOAD.replace("Info/Upload", chip)
                break
            else:
                PATH = ""

    return PATH


##############################     建立对应项目名的数据目录    ##################################
def mkdir(project):
    DATE = str(datetime.date.today()).replace("-", "")
    os.system("mkdir -p %s%s_%s_MGISEQ200" % (Sequences, project, DATE))
    target = "%s%s_%s_MGISEQ200" % (Sequences, project, DATE)
    return target


#######################################   数据拷贝   #########################################
def copy(PATH, target, project):
    os.system("cp -r %s %s" % (PATH, target))
    # os.chdir(TMP)
    # os.system("mv sample.list %s.xls" % project)
    with open(UnAnalysis,"a") as recored:
        recored.write("%s.xls\n" % project)
    return 1


#######################################     发送邮件    #########################################
def sendmail(done, target):
    if done:
        sender = 'kong.lulu@rigen-bio.com'
        receiver = ['kong.lulu@rigen-bio.com']
        message = MIMEText("数据路径：%s" % target, 'plain', 'utf-8')
        message['From'] = Header("研发项目", 'utf-8')  # 发送者
        message['To'] = Header("研发项目", 'utf-8')  # 接收者
        message['Subject'] = Header("【注册】下机数据已传输完成。", 'utf-8')  # 标题
        stp = smtplib.SMTP()
        stp.connect("smtp.mxhichina.com", 25)
        # stp.set_debuglevel(1)
        stp.login(sender, "822939Kll")
        stp.sendmail(sender, receiver, message.as_string())
        print("邮件发送成功")
        stp.quit()


##############################################################################################
def run():
    ready = diff(UnAnalysis,TotalTrans)
    for i in range(len(ready)):
        if os.path.exists(sample_dir+ready[i]):
            project, chip = load_sample(sample_dir+ready[i])
            if chip != "Not":
                PATH1 = Log_compare(UPLOAD1, chip)
                PATH2 = Log_compare(UPLOAD2, chip)
                PATH3 = Log_compare(UPLOAD3, chip)
                PATH4 = Log_compare(UPLOAD4, chip)
                PATH5 = Log_compare(UPLOAD5, chip)
                PATH = max(PATH1, PATH2, PATH3, PATH4, PATH5)
                if PATH == '':
                    print("项目:%s_芯片号:%s——数据未下机"%(project,chip))
                else:
                    target = mkdir(project)
                    done = copy(PATH, target, project)
                    sendmail(done, target)
                    break
            else:	#not zhuce xls
                with open(UnAnalysis,"a") as recored:
                    recored.write("%s.xls\n" % project)


if __name__ == '__main__':
    run()
