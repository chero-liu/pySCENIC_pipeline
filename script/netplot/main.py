# Author:liuchenglong

import argparse
import os
import sys

def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--regulons', help='rds file', required='T')
    parser.add_argument('-t', '--TF', help='', required='T')
    parser.add_argument('--highGene', help='', default='F')
    parser.add_argument('--legend', help='', default='T')

    parser.add_argument('--top', help='default 30', default=30)
    parser.add_argument('--showlabel', default=0.3, help='default 0.3')
    parser.add_argument('--vertexLabelSize', default=0.8, help='default 0.8')

    parser.add_argument('-o', '--outdir', help='outdir', default='.', required='T')
    parser.add_argument('-p', '--prefix', help='prefix', default='cellphone', required='T')

    parser.add_argument('--highGeneSize', help='', default=3)
    parser.add_argument('--highGeneColor', help='', default="00AA55")
    parser.add_argument('--TFColor', help='', default="1E90FF")
    parser.add_argument('--geneColor', help='', default="E87D72")
    parser.add_argument('--edgeColor', help='', default="00CED1")
    
    args = parser.parse_args()
    return(args)

def createDir(outdir, rm=False):
    if type(outdir).__name__ == 'list':
        for out in outdir:
            if not os.path.exists(out): 
                os.makedirs(out)
            else:
                if rm:
                    os.system(f"rm -r {out}")
                    os.makedirs(out)
    else:
        if not os.path.exists(outdir): 
            os.makedirs(outdir)
        else:
            if rm:
                os.system(f"rm -r {outdir}")
                os.makedirs(outdir)

if __name__ == '__main__':
    args = getargs()
    createDir(args.outdir, rm=False)
    string = f'''Rscript /SGRNJ03/PiplineTest01/Software_test/liuchenglong/pyscenicNetplot/utils/getdata.R --highGeneSize {args.highGeneSize} --geneColor {args.geneColor} --TFColor {args.TFColor} --highGeneColor {args.highGeneColor} --prefix {args.prefix} --outdir {args.outdir} --regulons {args.regulons} --TF {args.TF} --showlabel {args.showlabel} --top {args.top} --highGene {args.highGene}
Rscript /SGRNJ03/PiplineTest01/Software_test/liuchenglong/pyscenicNetplot/utils/netPlot.R --vertexLabelSize {args.vertexLabelSize} --edgeColor {args.edgeColor} --legend {args.legend} --geneColor {args.geneColor} --TFColor {args.TFColor} --highGeneColor {args.highGeneColor} --outdir {args.outdir} --prefix {args.prefix} --showlabel {args.showlabel} --highGene {args.highGene}
'''
    
    with open(f'{args.outdir}/pysenicNetPlot.sh', 'w') as sh:
        sh.write(string)
    os.system(string)
