# -*- coding: utf-8 -*-

import os
import re
import warnings
import unicodedata
import smtplib
from email.mime.text import MIMEText
from email.header import Header

warnings.filterwarnings("ignore")

'''
Usage: python main.py
'''
#######################    根据送样信息表的判断标准    ################################
#   建库试剂盒（判断是公司下机的还是明码的，下机数据格式不同）
#   建库试剂盒，样本类型，panel（1、A2，thyroscan，注册；2、DNA还是RNA；3、version）
#   下机数据目录（注册有"_MGISEQ200"后缀，没有后缀的话暂时建立A2的分析目录）
#	之后想加上：一个送样单里有几个项目，能建立对应的分析目录和分析脚本-----1
#		    分析脚本备份，能按不同的项目号和对应分析版本命名-----2
#################################   文件路径    #####################################
Analysis = "/data01/home/robot02/Project/2023/"  # 分析目录
Sequences = "/data01/project/sequence/2023_data/"  # 数据路径
sample_dir = "/data01/home/robot02/SampleList/"  # 送样信息表目录

Log = "/data01/home/robot02/.TMP/sequence.log"  # 备份下机数据log文件
sh_temp = "/data01/home/robot02/.TMP/_temp.sh"  # 分析脚本备份文件
Total = "/data01/home/robot02/.TMP/Total.log"  # 所有项目的送样信息表
Analysed = "/data01/home/robot02/.TMP/Analysed.log"  # 已提交分析的项目


###################################  查找R1 R2   ######################################
def find_file(filename, path):
    result = []
    for root, list, files in os.walk(path):
        for file in files:
            if filename in file:
                write = os.path.join(root, file)
                result.append(write)
    return result


#################################   下载数据 log    ################################
def log():
    os.system("ls %s > %s" % (Sequences, Log))
    os.system("ls %s > %s" % (sample_dir, Total))


###################################  找到没有分析的   ####################################
def diff(Total, Analysed):
    file_1 = open(Total, 'r')
    file_2 = open(Analysed, 'r')
    a = file_1.read().splitlines()
    b = file_2.read().splitlines()
    file_1.close()
    file_2.close()
    ready = ([x for x in a if x not in b])
    if len(ready) == 0:
        raise SystemExit('没有待分析!')
    else:
        return ready


###############################     分析脚本备份文件    ##############################
def runtmp(samplelist):
    dir = {}
    with open(sh_temp, 'w') as ftmp:
        ################################     读送样信息表    ################################
        with open(samplelist, 'r') as f:
            for line in f:
                r1 = ""
                r2 = ""
                line = line.strip('\n').split('\t')
                samplename = line[0] + "-" + line[1]
                ProjectNumber = line[0]
                type = line[4]
                sample = line[1]
                version = unicodedata.normalize('NFKC', line[6])
                ##########################      对应项目的数据路径     ##############################
                path = ""
                with open(Log, "r") as flog:
                    for log in flog:
                        log = log.strip('\n')
                        searchObj = re.search(ProjectNumber, log)
                        if searchObj:
                            path = Sequences + log
                ##########################    根据样本名找到对应的R1,R2   ###########################
                ####################   明码下机  &&  公司内测序   #####################
                if "RJ-YZ" in line[0]:  # 公司的
                    total = find_file(sample, path)
                    for i in range(len(total)):
                        if "L01_"+sample+"_1.fq.gz" in total[i]:
                            r1 = total[i]
                        if "L01_"+sample+"_2.fq.gz" in total[i]:
                            r2 = total[i]
                else:  # 明码的
                    total = find_file(samplename, path)
                    for i in range(len(total)):
                        if samplename + "_combined_R1.fastq.gz" in total[i] and "md5" not in total[i]:
                            r1 = total[i]
                        if samplename + "_combined_R2.fastq.gz" in total[i] and "md5" not in total[i]:
                            r2 = total[i]

                ################################   项目及版本  ##################################

                A2 = {"DNA": {
                    "3.3": "perl /data01/home/robot02/thyroid_18gene_A2/DNA/Main_program/main_amp_dna_A2_v3.3.pl -1 {0:s} -2 {1:s} -s {2:s} -p ./".format(
                        r1, r2, samplename)},
                    "RNA": {
                        "1.2": "perl /data01/home/robot02/thyroid_18gene_A2/RNA/Main_program/main_amp_rna_A2_v1.2.pl -1 {0:s} -2 {1:s} -s {2:s} -p ./".format(
                            r1, r2, samplename),
                        "1.3": "perl /data01/home/robot02/thyroid_18gene_A2/RNA/Main_program/main_amp_rna_A2_v1.3.pl -1 {0:s} -2 {1:s} -s {2:s} -p ./".format(
                            r1, r2, samplename)
                    }}


                zhuce = {"DNA": {
                    "1.1": "perl /data01/home/robot02/zhuce_6gene/DNA/Main_program/main_v2.pl -1 {0:s} -2 {1:s} -s {2:s} -p ./ -c /data01/home/robot02/zhuce_6gene/DNA/config/thyroid_ZC_v4_short.cfg -r".format(
                        r1, r2, samplename),
                    "1.3": "perl /data01/home/robot02/zhuce_6gene/DNA/Main_program/main_v2.pl -1 {0:s} -2 {1:s} -s {2:s} -p ./ -c /data01/home/robot02/zhuce_6gene/DNA/config/thyroid_ZC_v2.0_short.cfg -r".format(
                        r1, r2, samplename),
                    "2.1": "perl /data01/home/robot02/zhuce_6gene/DNA/Main_program/main_v2.pl -1 {0:s} -2 {1:s} -s {2:s} -p ./ -c /data01/home/robot02/zhuce_6gene/DNA/config/thyroid_ZC_v2.1.cfg -r".format(
                        r1, r2, samplename),
                    "2.2": "perl /data01/home/robot02/zhuce_6gene/DNA/Main_program/main_v2.pl -1 {0:s} -2 {1:s} -s {2:s} -p ./ -c /data01/home/robot02/zhuce_6gene/DNA/config/thyroid_ZC_v2.2.cfg -r".format(
                        r1, r2, samplename), },
                    "RNA": {
                        "1.0": "perl /data01/home/robot02/zhuce_6gene/RNA/Main_program/main_rna_v2.pl -1 {0:s} -2 {1:s} -s {2:s} -p ./ -c /data01/home/robot02/zhuce_6gene/RNA/config/amp_rna_ZC_v4.cfg".format(
                            r1, r2, samplename)}}
                ###############################    匹配上对应主程序    ##############################
                my_command = ""
                ver = ''
                if "RJ-YA" in line[0]:
                    for i in A2[type].keys():
                        if i in version:
                            ver = i
                            my_command = A2[type][ver]
                elif "RJ-YZ" in line[0]:
                    for i in zhuce[type].keys():
                        if i in version:
                            ver = i
                            my_command = zhuce[type][ver]
                else:
                    print("ERROR: 项目 %s 的样本类型/版本号/项目类型可能有误。" % version)
                    print(line)
                dir[type] = ver
                ftmp.write("%s\n" % (my_command))
    return ProjectNumber, dir, path


##################################    建立分析目录   ##################################
def mkdir_anly(dir, path):
    if "_MGISEQ200" in path.split("/")[5]:
        if len(dir) == 2:
            PATH = path.split("/")[5].replace("_MGISEQ200",
                                              "_DNA_Hotspot_v%s_RNA_fusion_v%s" % (dir["DNA"], dir["RNA"]))
        elif len(dir) == 1:
            if "DNA" in dir.keys():
                PATH = path.split("/")[5].replace("_MGISEQ200", "_DNA_Hotspot_v%s" % dir["DNA"])
            else:
                PATH = path.split("/")[5].replace("_MGISEQ200", "_RNA_fusion_v%s" % dir["RNA"])
        else:
            print("样本类型错误！")
    else:
        if len(dir) == 2:
            PATH = path.split("/")[5] + "_A2_DNA_v%s_RNA_v%s" % (dir["DNA"], dir["RNA"])
        elif len(dir) == 1:
            if "DNA" in dir.keys():
                PATH = path.split("/")[5] + "_A2_DNA_v%s" % dir["DNA"]
            else:
                PATH = path.split("/")[5] + "_A2_RNA_v%s" % dir["RNA"]
        else:
            print("样本类型错误！")
    Path = Analysis + PATH
    os.system("mkdir -p %s" % Path)
    return Path


##################################    运行    ####################################
def copy_run(Path, ready):
    my_command1 = "cp %s %s/run.sh" % (sh_temp, Path)
    os.system(my_command1)
    with open(Analysed, "a") as fan:
        fan.write("%s\n" % ready)
    os.chdir(Path)
    os.system("sh ./run.sh")
    return 1


####################################发邮件##########################################
def sendmail(done, path, ProjectNumber, Path):
    if done:
        sender = 'kong.lulu@rigen-bio.com'
        receiver = ['kong.lulu@rigen-bio.com','qiao.congcong@rigen-bio.com']
        message = MIMEText("分析路径：%s" % Path, 'plain', 'utf-8')
        message['From'] = Header("研发<Analy>", 'utf-8')  # 发送者
        message['To'] = Header("研发<Analy>", 'utf-8')  # 接收者
        if "_MGISEQ200" in path.split("/")[5]:
            message['Subject'] = Header("【注册项目】 %s 已提交分析。" % ProjectNumber, 'utf-8')  # 标题
        else:
            message['Subject'] = Header("【多重18基因A2项目】 %s 已提交分析。" % ProjectNumber, 'utf-8')  # 标题
        stp = smtplib.SMTP()
        stp.connect("smtp.mxhichina.com", 25)
        # stp.set_debuglevel(1)
        stp.login(sender, "822939Kll")
        stp.sendmail(sender, receiver, message.as_string())
        print("邮件发送成功")
        stp.quit()


##################################################################################
def run():
    log()
    ready = diff(Total, Analysed)
    for i in range(len(ready)):
        ProjectNumber, dir, path = runtmp(sample_dir + ready[i])
        if path != "":
            Path = mkdir_anly(dir, path)
            done = copy_run(Path, ready[i])
            sendmail(done, path, ProjectNumber, Path)
            break
        else:
            print("项目:%s——数据未下机。"%ProjectNumber)

if __name__ == '__main__':
    run()
