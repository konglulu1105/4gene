import os
import smtplib
from email.mime.text import MIMEText
from email.header import Header
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart

Analysis = "/data01/home/robot02/Project/2023/"  # 分析目录
sample_dir = "/data01/home/robot02/SampleList/"  # 送样信息表目录
python = "/data01/software/miniconda3/envs/py3/bin/python"
Script = "/data01/home/robot02/script/Main_report.py"
Analysed = "/data01/home/robot02/.TMP/Analysed.log"  # 提交了分析的项目
Reported = "/data01/home/robot02/.TMP/Reported.log"  # 汇总了的分析结果


########################################读取已分析的项目表###########################################
def diff(Analysed, Reported):
    file_1 = open(Analysed, 'r')
    file_2 = open(Reported, 'r')
    a = file_1.read().splitlines()
    b = file_2.read().splitlines()
    file_1.close()
    file_2.close()
    ready = ([x for x in a if x not in b])
    if len(ready) == 0:
        raise SystemExit('没有待报告!')
    else:
        return ready

###########################################样本信息表数量############################################
def Samp_num(sample_list):
    hang = len(open(sample_list).readlines())
    return hang


###########################################    QC数量    ########################################
def QC_num(path):
    file = os.listdir(os.path.join(path,"results"))
    qc_num = sum([target.count("QC.xls") for target in file])
    return qc_num


#############################################找到对应分析目录########################################
def find_file(ready, Analysis):
    files = os.listdir(Analysis)
    #print(Analysis)
    target = ready.removesuffix(".xls")
    for i in files:
    #    print(i) 
        if target in i:
            path = Analysis + i
    #print(path)
    return path


##############################################汇总#############################################
def make_rep(project, path):
    os.chdir(path)
    proj = project.removesuffix(".xls")
    meslog = os.popen("%s %s %s" % (python, Script, proj))
    with open(Reported,"a") as frep:
        frep.write("%s\n"%project)
    files = os.listdir(path)
    file = []
    for i in files:
        if i.endswith("xls"):
            file.append(i)
    return 1, file, meslog


#############################################sendemail##########################################
def sendmail(done, path, files, project, meslog):
    if done:
        sender = 'kong.lulu@rigen-bio.com'
        receiver = ['kong.lulu@rigen-bio.com','qiao.congcong@rigen-bio.com','huang.yan@rigen-bio.com']
        message = MIMEMultipart()
        part = MIMEText("分析路径：%s\n%s" % (path,meslog.read()), 'plain', 'utf-8')
        message.attach(part)
        message['From'] = Header("研发项目", 'utf-8')  # 发送者
        message['To'] = Header("研发项目", 'utf-8')  # 接收者
        if "A2" in path:
            message['Subject'] = Header("【多重18基因A2项目】 %s 已分析完成。" % project, 'utf-8')  # 标题
        else:
            message['Subject'] = Header("【注册项目】 %s 已分析完成。" % project, 'utf-8')  # 标题

        for i in files:
            part = MIMEApplication(open(i, 'rb').read())
            name = i.replace(path, "")
            part.add_header('Content-Disposition', 'attachment', filename=name)
            message.attach(part)
        stp = smtplib.SMTP()
        stp.connect("smtp.mxhichina.com", 25)
        # stp.set_debuglevel(1)
        stp.login(sender, "822939Kll")
        stp.sendmail(sender, receiver, message.as_string())
        print("邮件发送成功")
        stp.quit()


def run():
    ready = diff(Analysed,Reported)
    for i in range(len(ready)):
        samp_num = Samp_num(os.path.join(sample_dir,ready[i]))
        print(ready[i])
        path = find_file(ready[i], Analysis)
        print(path)
        qc_num = QC_num(path)
        if samp_num == qc_num:
            done, file, meslog = make_rep(ready[i], path)
            sendmail(done, path, file, ready[i].removesuffix(".xls"), meslog)
            break
if __name__ == '__main__':
    run()
