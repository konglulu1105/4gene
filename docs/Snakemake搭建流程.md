# Snakemake搭建流程
## 介绍

对于生信分析, 大部分人的第一印象就是流程化。确实，再依托于强大的计算机系统，许多常规的分析任务都可以自动化运行并进行相应的分析。例如对于测序数据来说， 常见的上游分析便是数据质控、去接头、序列比对、PCR去重矫正、基因表达定量、突变和拷贝数分析等。这些分析步骤中，每一步都有许多不同的软件和算法，选择多样化。而对于数据分析来说，这些上游的步骤并不是他们所关注的，他们可能更注重一些下游的个性化分析，所以前面的数据处理和分析工作完全可以让流程自动完成。

所谓流程，便是能够将分析步骤进行串联，多样本并行计算。常见的流程化一般使用`shell`脚本，将所有分析任务写入到一个脚本中，这种方法扩展性较差，程序臃肿。而`Python`是一门非常灵活的语言，非常适合用来编写模块化的流程，但是比较麻烦的一点是，需要自己进行一些命令的封装。幸好，已经有人帮我们做好了。

`Snakemake`是一种基于`Python`的工作流管理系统，能够用于创建可重复和可扩展的数据分析的工具。其编写的流程可以无缝部署到服务器、集群等环境，而无需对流程进行修改。

## 定义流程

受`Make`的启发，`snakemake`使用`rule`字段来定义一个从输入文件到输出文件的规则，而不同规则之间的依赖关系是通过输入和输出文件名称继续匹配，隐式进行处理。一个`rule`可以看作是一个操作步骤，必须序列比对、表达定量等调用特定软件的步骤就可以定义一个`rule`来处理。

### 1. 基本语法

-   `snakemake`使用冒号和缩进来定义规则，并用一些关键字来进行规则声明，同时也可以兼容`Python`代码
-   在`rule`声明的规则中，有很多关键字用于定义不同的规则和文件信息，例如

| 关键字 | 作用                         | 关键字 | 作用         |
| ------ | ---------------------------- | ------ | ------------ |
| input  | 定义输入文件                 | output | 定义输出文件 |
| params | 定义执行命令时用到的额外参数 | log    | 定义日志文件 |
| threads |	定义该规则使用的线程数	|	message	|打印消息|
|resouces|	定义使用资源，如内存，CPU等	| version |	定义规则的版本	|
| conda | 定义conda环境 | container | 运行在容器环境中 |
| run | 定义多行命令 | shell | 执行shell命令 |
| script | 定义规则中的处理脚本 | notebook | 执行jupyter notebook文件 | 

### 2. 定义`rule`

- 一个`rule`定义流程中的一项数据分析步骤，使用关键字`rule`声明，包括规则名称，输入文件、输出文件以及输入文件映射到输出文件的`shell`命令或一个处理脚本。例如

```python
rule Name:
	input: "path/to/inputfile1", "path/to/inputfile2"
	output: ["path/to/outputfile2", "path/to/outputfile2"]
	shell: "command {input} {output}"
```

-   规则的名称是可选的，可以省略，表示一个匿名规则。也可以通过设置规则的`name`属性来覆盖。
-   其中`input`和`output`可以是元组或者列表，最后一行`shell`表示需要执行的`shell`命令，其中花括号表示引用，即使用`input`和`output`关键字中定义的值来进行字符串替换。可以使用双花括号来取消这种占位符替换，会被转换为单花括号字符。在这里输入和输出列表会被转换为空格分隔的字符串（类似于 `''.join(input)`）
-   也可以使用字典的方式定义键值对，例如

```python
rule NAME:
    input: 
        fastq1 = "path/to/inputfile", 
        fastq2 = "path/to/other/inputfile"
    output: "path/to/outputfile", somename = "path/to/another/outputfile"
    shell: "command {input.fastq1} {input.fastq2} {output[0]}"
```

-   在花括号内使用点加属性名称或使用数组索引的方式访问
-   `TODO` `shell`命令描述



## 通配符

1.   有时候， 我们的输入或输出文件并不是固定的名称，比如我们分析一批数据不可能只有一个样本，而且每次分析样本的名称肯定不一样，因此将输入或输出指定为一种可变的，或可替换的字符是很有必要的，这也是`snakemake`的一个优势所在。例如下面这个例子。

```python
rule complex_conversion:
    input:
        "input/{sample_name}_inputfile.txt"
    output:
        "output/{sample_name}_outputfile.txt"
    shell:
        "somecommand --sample_name {wildcards.sample_name} < {input} > {output}"
```

-   在这里， 我们定义一个通配符`sample_name`，这个通配符相当于正则表达式中的`.+`，可以匹配除空字符外的任意字符串。

-   在`shell`字段中，可以使用`wildcards`对象来访问对应的通配符，该对象是个内置对象，所有通配符均为其属性。

2.   `rule`规则定义在`Snakefile`中。也可以在`rule all`中`input`定义。

     ```python
     rule all:
         input:
             "output/A/file.R1.txt"
     ```

     -   `all`这个规则名称是固定的，而且只要在`input`中定义执行的输出文件，`snakemake`会自动运行相应的规则，及其依赖的规则。也就是说如果该文件的输出需要其他规则，那些规则也会自动执行，其最终的目的也就必须获取到`input`里面定义的文件
     -   `snakemake`默认指挥

