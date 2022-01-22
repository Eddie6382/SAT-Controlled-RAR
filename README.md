# SAT-Controlled RAR

## reference
1. C. -A. Wu, T. -H. Lin, Shao-Lun Huang and C. -Y. Huang, "SAT-Controlled redundancy addition and removal — a novel circuit restructuring technique," 2009 Asia and South Pacific Design Automation Conference, 2009, pp. 191-196, doi: 10.1109/ASPDAC.2009.4796479.

2. DSNP course by Prof. Ric

## requirement
python3 for ```checker.py```  
c++ > 11 for ```cirTest```  

## Usage
./cirTest  
```
   make clean   # clean
   make
   ./cirTest [-F] [doFile]
```

check.py  
```
   python3 check.py [args....]
```
```--check```: To verify each repairment or not, if a repaired circuit behaved not identical to the original circuit, it will existed in directory ```results/err_report```  
```--benchmark_dir```  
```--file_name```: specify one file in benchmark_dir  
```--result_dir```: dump logout to result_dir  

## Examples
**mydo**  
```
   cirr ISCAS85/ex01.aag
   cirrar -v 1
   cirrarw -directory results
   q -f
```
**results**
```
 /mnt/e/u/O/n/Ric_SAT/SATRar | master *1 !12  ./cirTest -F mydo                                                ok | 01:14:02
cir> cirr ISCAS85/ex01.aag

cir> cirrar -verb 0

----  Preprocess ----
Done
----  Run SATrar ----
  w_t(8, 11)
  w_t(11, 14)
  w_t(9, 12)
  w_t(12, 13)
        conflict gid 9
        conflict gid 9
  w_t(10, 13)
  w_t(13, 14)
  w_t(14, 7)

----  Statistic  ----
Total time usage: 0 s
#wire: 7
2-Way rar alternatives: 0
RAR wire alternatives: 0
RAR gate alternatives: 2
#tar: 1
#alt: 2

cir> cirrarw -directory results
writing files...done

cir> q -f
```

in ```results/repair0.aag```
```
aag 17 6 0 1 10
2
4
6
8
10
12
14
16 4 8
22 6 16
18 10 7
24 19 9
20 2 4
32 1 20
34 17 19
26 32 35
28 23 27
14 29 12
c
aag repair, redundant wire = (12, 13)
            g_d = 13
            conflict decisions = 8 9 
```

**vertification**  
```
/mnt/e/u/O/n/Ric_SAT/SATRar | master *1 !15  python3 check.py --check --benchmark_dir ISCAS85 --file_name ex01.aag --result_dir results
SUCCESS: ISCAS85/ex01.aag
    Correctness check on repair0.aag... success
    Correctness check on repair1.aag... success
```

