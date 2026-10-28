"""Build the optical report from archived results; never rerun either model.

python scripts/build_optical_validation_report.py --charts
python scripts/build_optical_validation_report.py --pdf
The first step needs matplotlib; the second needs reportlab, Pillow and pypdf.
"""
import argparse
import csv
import hashlib
import json
import math
import shutil
from pathlib import Path
from xml.sax.saxutils import escape

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'docs/validation/report_20261024'
BUNDLE = ROOT / 'docs/validation/current_model_20261021'
OUT = ROOT / 'docs/validation/optical_report_20261024'
MD = ROOT / 'docs/TIRT_光学验证报告.md'
PDF = ROOT / 'docs/TIRT_光学验证报告.pdf'
LABELS = ['平地裸土', '平地均质冠层', '平地离散森林', '垄行（垄向 0°）', '垄行（垄向 90°）',
          '单坡＋均质冠层', '单坡＋森林', '复合山地＋森林', '城市（无植被）', '城市＋森林']


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    with path.open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stats(rows):
    differences = [float(r['value']) - float(r['reference']) for r in rows]
    assert differences and all(math.isfinite(v) for v in differences)
    return dict(n=len(rows), rmse=math.sqrt(sum(d*d for d in differences)/len(rows)),
                bias=sum(differences)/len(rows), max_abs=max(map(abs, differences)))


def fmt(value, signed=False):
    value = float(value)
    if value and abs(value) < 1e-4:
        return f'{value:+.2e}' if signed else f'{value:.2e}'
    return f'{value:+.4f}' if signed else f'{value:.4f}'


def load():
    sources = ['accuracy.csv', 'reference_comparison.csv', 'verification.json',
               'scenario_verification.json', 'optical_lai_spectrum.csv', 'optical_lai_angles.csv',
               'optical_lai_spectrum.png', 'optical_lai_angles.png', 'optical_rmse.png']
    hashes = read_json(SOURCE/'report_manifest.json')['files_sha256']
    for name in sources:
        path = SOURCE/name
        assert sha(path) == hashes[str(path.relative_to(ROOT)).replace('/', '\\')], name
    scenes = read_json(BUNDLE/'manifest.json')['scenes']
    paired = [r for r in read_csv(SOURCE/'reference_comparison.csv') if int(r['wavelength_nm']) in (670, 860)]
    accuracy = [r for r in read_csv(SOURCE/'accuracy.csv') if int(r['wavelength_nm']) in (670, 860)]
    assert len(paired) == 1204 and len(scenes) == 10 and len(accuracy) == 20
    assert len({(r['scene'], r['wavelength_nm'], r['vza'], r['vaa']) for r in paired}) == 1204
    scores = {}
    for row in accuracy:
        subset = [r for r in paired if r['scene'] == row['scene'] and r['wavelength_nm'] == row['wavelength_nm']]
        score = stats(subset)
        for key, value in score.items():
            assert math.isclose(value, float(row[key]), abs_tol=1e-12, rel_tol=1e-12), (row['scene'], key)
        assert score['n'] == (57 if 'slope' in row['scene'] else 61)
        scores[row['scene'], int(row['wavelength_nm'])] = score
    for row in paired:
        assert abs(float(row['difference']) - (float(row['value'])-float(row['reference']))) < 1e-12
    OUT.mkdir(parents=True, exist_ok=True)
    write_csv(OUT/'reference_comparison.csv', paired)
    write_csv(OUT/'accuracy.csv', accuracy)
    totals = [dict(wavelength_nm=w, **stats([r for r in paired if int(r['wavelength_nm']) == w])) for w in (670, 860)]
    write_csv(OUT/'overall_accuracy.csv', totals)
    coverage = [r for r in read_json(SOURCE/'scenario_verification.json')['details'] if r['domain'] == 'optical']
    assert len(coverage) == 10 and sum(r['valid'] for r in coverage) == 858
    write_csv(OUT/'scenario_coverage.csv', coverage)
    for name in sources[4:]:
        shutil.copy2(SOURCE/name, OUT/name)
    source_paths = [SOURCE/name for name in sources] + [BUNDLE/'manifest.json', BUNDLE/'reference_manifest.json',
                    BUNDLE/'histream_reference.csv'] + [BUNDLE/'inputs'/f"{s['stem']}.json" for s in scenes]
    (OUT/'source_manifest.json').write_text(json.dumps({
        'report_date': '2026-10-24', 'research_source': '闻建光团队',
        'result_snapshot': 'report_20261024', 'input_bundle': 'current_model_20261021',
        'report_build_reran_tirt': False, 'report_build_reran_histream': False,
        'scene_count': 10, 'optical_pairs': 1204, 'metrics_recalculated_from_pairs': True,
        'files': {p.relative_to(ROOT).as_posix(): sha(p) for p in source_paths}
    }, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    return scenes, paired, scores, totals, coverage


def charts(scenes, paired):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.sans-serif': ['Microsoft YaHei', 'DejaVu Sans'],
                         'axes.unicode_minus': False, 'font.size': 10})
    for scene, label in zip(scenes, LABELS):
        fig, panels = plt.subplots(2, 2, figsize=(10, 7.2), sharex=True, sharey='row')
        for row_index, wave in enumerate((670, 860)):
            rows = [r for r in paired if r['scene'] == scene['stem'] and int(r['wavelength_nm']) == wave]
            for col, (positive, negative) in enumerate(((0, 180), (90, 270))):
                ax = panels[row_index, col]
                selected = sorted([(0. if float(r['vza']) == 0 else float(r['vza'])*(1 if float(r['vaa']) == positive else -1), r)
                                  for r in rows if float(r['vza']) == 0 or float(r['vaa']) in (positive, negative)], key=lambda v:v[0])
                for key, name, color, marker in [('reference', 'HiSTREAM', '#e68a2e', 'o'), ('value', 'TIRT', '#167c89', None)]:
                    ax.plot([v[0] for v in selected], [float(v[1][key]) for v in selected],
                            label=name, color=color, marker=marker, markersize=3, linewidth=1.6)
                ax.set_title(f"{wave} nm | {'太阳主平面' if col == 0 else '垂直太阳主平面'}")
                ax.set_xlim(-80, 80)
                ax.grid(alpha=.22)
                ax.set_ylabel('方向反射因子 BRF' if col == 0 else '')
                ax.set_xlabel('带符号观测天顶角 (°)')
                ax.legend(loc='best')
            values = [float(r[k]) for r in rows for k in ('value', 'reference')]
            pad = max(.07*(max(values)-min(values)), .005)
            panels[row_index, 0].set_ylim(max(0, min(values)-pad), max(values)+pad)
        fig.suptitle(label+' | SZA 30°, SAA 0°', fontsize=14)
        fig.tight_layout()
        fig.savefig(OUT/f"{scene['stem']}.png", dpi=190)
        plt.close(fig)


def document(scenes, scores, totals, coverage):
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak
    from PIL import Image as PILImage
    from pypdf import PdfReader
    pages = []
    def p(value): return ('text', value)
    def sub(value): return ('subheading', value)
    def table(header, rows, widths): return ('table', (header, rows, widths))
    def figure(name, caption, height): return ('figure', (name, caption, height))
    def page(title, *blocks): pages.append([('heading', title), *blocks])
    page('TIRT 光学验证报告',
         p('报告日期：2026-10-24 | 10 个场景 | 红光与近红外方向反射因子'),
         sub('研究来源与报告依据'),
         p('TIRT 光学部分来自闻建光团队的相关研究工作。本报告依据仓库已保存的 2026-10-24 TIRT 复算结果和冻结 HiSTREAM 参考结果整理，不是对原始论文结果的逐项复现。'),
         p('本次报告编制重新核算逐方向误差、绘制曲线和整理运行示例，未重新运行 TIRT 或 HiSTREAM。输入场景沿用 current_model_20261021 数据包。与热红外主页报告的版本和场景统计口径分别记录。'),
         sub('验证范围与总体结果'),
         p('覆盖平地裸土、均质冠层、离散森林、两种垄向作物、单坡植被、复合山地森林及城市植被。每波段 602 条有效配对，共 1,204 条；相对 HiSTREAM 的误差单位为无量纲 BRF。'),
         table(['波段', 'n', '整体 RMSE', '整体偏差', '最大绝对误差'],
               [[str(r['wavelength_nm'])+' nm', str(r['n']), *[f"{r[k]:.6f}" for k in ('rmse', 'bias', 'max_abs')]] for r in totals], [100,45,110,110,140]),
         p('偏差 = TIRT - HiSTREAM；RMSE = 差值平方均值的平方根。整体指标按全部有效配对计算，不是场景 RMSE 的简单平均。'),
         sub('照明、观测和材料设置'),
         p('太阳天顶角 30°、方位角 0°，太阳直射，漫射比例为 0。观测天顶角 0° 至 75°、步长 5°；方位角为 0°、90°、180°、270°，天底只统计一次。通常每场景每波段 61 个方向，两个单坡场景排除背坡后各为 57 个。'),
         table(['材料参数', '670 nm', '860 nm'], [['叶片反射率','0.06','0.45'],['叶片透射率','0.03','0.42'],['土壤 / 街道反射率','0.18','0.30'],['屋顶反射率','0.22','0.35'],['墙面反射率','0.30','0.45']], [250,125,125]),
         p('不运行能量平衡。BRF 是方向反射因子，不是半球反照率；本组结果是模型间对比，不能等同于地面或卫星观测精度。'))
    rows = []
    for scene, label in zip(scenes, LABELS):
        red, nir = (scores[scene['stem'], w] for w in (670,860))
        rows.append([label,str(red['n']), *[fmt(s[k], k=='bias') for s in (red,nir) for k in ('rmse','bias','max_abs')]])
    page('1　分场景精度汇总',
         p('研究基础：闻建光团队。n 为每波段有效方向数；所有误差均相对冻结 HiSTREAM，单位 BRF。'),
         table(['场景','n','红光\nRMSE','红光\n偏差','红光\n最大误差','近红外\nRMSE','近红外\n偏差','近红外\n最大误差'], rows,[120,27,65,68,70,65,68,70]),
         figure('optical_rmse.png','图 1　十个场景的红光与近红外 RMSE。',255),
         p('红光 RMSE 最大为城市＋森林（0.0098）；近红外 RMSE 最大为垄向 90° 作物（约 0.0217）。近红外最大逐方向误差为单坡＋森林的 0.084806，平均误差不能代替局部方向差异。'),
         p('裸地误差在 1e-7 量级，保留科学计数法，避免将舍入后的零解释为严格一致。'))
    notes = [
        '均匀裸土基准；两波段土壤反射率分别为 0.18、0.30，用于检查无植被条件下的角度响应。',
        '均质冠层 LAI 为 1.5。红光以吸收为主，近红外具有较强的叶片散射和透射。',
        '离散森林采用 0.15 m 叶片尺度的参考场景；冠层 LAI 为 1.5。',
        '垄宽 0.5 m、垄间空白 0.5 m、垄高 1 m，垄向为 0°。',
        '垄宽、垄间空白和垄高保持不变，垄向改为 90°；需同时查看两个观测主平面。',
        '30° 朝东单坡，上覆均质冠层。专用验证脚本先转换太阳和观测方向到坡面局地坐标，再调用平地光学模型。',
        '30° 朝东单坡，上覆离散森林；同样采用局地角度转换，背坡观测不进入误差统计。',
        '复合山地与森林耦合场景，山高与底面半径均为 100 m。验证脚本调用底层 directional_reflectance；与热红外主页中的 20 m 山体不同。',
        '建筑尺寸为 10 × 10 × 20 m，密度每平方米 0.0015 栋；街道不含植被。',
        '城市建筑与街道森林耦合，叶片尺度 0.15 m；光学输入保留植被、街道、屋顶和墙面的独立材料参数。']
    for index, (scene, label, note) in enumerate(zip(scenes,LABELS,notes),2):
        table_rows=[]
        for wave in (670,860):
            s=scores[scene['stem'],wave]
            table_rows.append([str(wave)+' nm',str(s['n']),fmt(s['rmse']),fmt(s['bias'],True),fmt(s['max_abs'])])
        page(f'{index}　{label}', p(note),
             table(['波段','n','RMSE','偏差','最大绝对误差'],table_rows,[95,40,115,115,140]),
             figure(scene['stem']+'.png',f'图 {index}　{label}：红光、近红外在两个主平面的 TIRT / HiSTREAM 对比。',420),
             p('每一波段的两个主平面采用相同纵轴范围。太阳主平面正角对应方位角 0°，负角对应 180°；垂直主平面正角对应 90°，负角对应 270°。天底在两幅图中均显示，统计时仅计一次。'),
             p('此曲线对比针对所列输入和方向。未设置统一的物理精度通过阈值；1e-6 仅用于数值复算的一致性检查。'))
    page('12　多 LAI、多光谱与多角度示例',
         p('平地均质植被，LAI = 1、3、5，太阳天顶角 30°。使用 optical_run.py 的已保存运行输出；材料光谱为案例设定值。'),
         figure('optical_lai_spectrum.png','图 12　天底方向下的七波段 BRF，横坐标为波长。',270),
         figure('optical_lai_angles.png','图 13　860 nm 的双主平面多角度 BRF，横坐标为观测天顶角。',220),
         p('七个示例波段为 450、550、670、705、860、1610、2200 nm；角度扫描为 -60° 至 60°。这些示例检验参数扫描与输出，只有 670、860 nm 另有本报告列出的参考模型对比。'))
    scene_names = dict(zip(['plane_bare','plane_hom','plane_row','plane_crown','slope_crown','terrain_crown','urban_bare','urban_hom','urban_row','urban_crown'],
                          ['平地裸土','平地均质','平地垄行','平地森林','单坡森林','山地森林','城市无植被','城市均质植被','城市垄行','城市森林']))
    page('13　运行覆盖、复现与适用范围',
         p('已有公开光学 JSON 入口测试共 10 个案例、858 条有效记录、20 张 PNG；每例保存 spectrum.csv、angles.csv 及两张图。以下为运行和输出完整性检查，不表示独立观测验证通过。'),
         table(['场景','有效记录','PNG','输出检查'],[[scene_names[r['scene']],str(r['valid']),str(r['pngs']),'通过'] for r in coverage],[245,95,65,100]),
         p('公开光学入口现已支持单坡 / 复合地形与裸土、均质、垄行的组合：坡面使用局地平面变换，复合地形使用组分比例。这些近似没有独立验证；单坡均质参考对比仍使用专用验证适配器。'),
         sub('证据与数值复现'),
         p('独立报告数据位于 docs/validation/optical_report_20261024/，包含 1,204 条逐方向配对、20 行分场景统计、双波段总体统计、图件及 LAI 示例 CSV。source_manifest.json 记录源文件 SHA-256；编制时从配对数据重新核算 RMSE、偏差与最大误差。'),
         p('原复算回执位于 report_20261024/verification.json：三波段联合复算采用绝对容差 1e-6、相对容差 0，最大差异 7.555229331e-9。该最大值包含热红外，本报告不将其误标为光学专属误差。'),
         p('参考来源由 current_model_20261021/reference_manifest.json 记录，其中光学源文件为 github_optical_reference.csv。现有 reference_protocols.json 主要记录 10.5 μm 参数，不能据此宣称光学也采用相同采样数与追踪深度。本次没有新增光线追踪采样不确定性评估。'),
         ('code','python validation/validate_current_model.py --plots\npython scripts/build_optical_validation_report.py --charts\npython scripts/build_optical_validation_report.py --pdf'),
         p('结论：本组光学结果覆盖 10 个典型结构场景的红光与近红外方向响应；城市森林、垄行及单坡森林仍需关注局部方向误差。研究来源：闻建光团队。'))
    text=[]
    for index, blocks in enumerate(pages):
        for kind, value in blocks:
            if kind=='heading': text.append(('# ' if index==0 else '## ')+value)
            elif kind=='subheading': text.append('### '+value)
            elif kind=='text': text.append(value)
            elif kind=='code': text.append('```powershell\n'+value+'\n```')
            elif kind=='table':
                header,values,_=value
                text.append('\n'.join(['| '+' | '.join(h.replace('\n',' ') for h in header)+' |','| '+' | '.join('---' for _ in header)+' |',*['| '+' | '.join(row)+' |' for row in values]]))
            elif kind=='figure':
                filename,caption,_=value
                text.append(f'![{caption}](validation/optical_report_20261024/{filename})\n\n{caption}')
    text.append('数据：[逐方向对比](validation/optical_report_20261024/reference_comparison.csv) · [精度统计](validation/optical_report_20261024/accuracy.csv) · [来源校验清单](validation/optical_report_20261024/source_manifest.json)')
    MD.write_text('\n\n'.join(text)+'\n',encoding='utf-8')
    pdfmetrics.registerFont(TTFont('CJK',r'C:\Windows\Fonts\simsun.ttc',subfontIndex=0))
    width=A4[0]-80
    styles={
        'heading':ParagraphStyle('h',fontName='CJK',fontSize=18,leading=25,spaceAfter=16,textColor=colors.HexColor('#153e4b')),
        'subheading':ParagraphStyle('s',fontName='CJK',fontSize=12,leading=17,spaceBefore=8,spaceAfter=6,textColor=colors.HexColor('#167c89')),
        'text':ParagraphStyle('b',fontName='CJK',fontSize=10,leading=16,spaceAfter=10,wordWrap='CJK'),
        'cell':ParagraphStyle('c',fontName='CJK',fontSize=8.5,leading=12,wordWrap='CJK'),
        'caption':ParagraphStyle('cap',fontName='CJK',fontSize=9,leading=13,spaceAfter=10,textColor=colors.HexColor('#50626b')),
        'code':ParagraphStyle('code',fontName='Courier',fontSize=8,leading=12,spaceAfter=10)}
    story=[]
    for index,blocks in enumerate(pages):
        if index: story.append(PageBreak())
        for kind,value in blocks:
            if kind in styles:
                story.append(Paragraph(escape(value).replace('\n','<br/>'),styles[kind]))
            elif kind=='table':
                header,values,widths=value
                data=[[Paragraph(escape(v).replace('\n','<br/>'),styles['cell']) for v in row] for row in [header,*values]]
                t=Table(data,colWidths=[w*width/sum(widths) for w in widths],repeatRows=1)
                t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#dfedf0')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f5f8f9')]),('VALIGN',(0,0),(-1,-1),'MIDDLE'),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]))
                story.extend([t,Spacer(1,10)])
            elif kind=='figure':
                filename,caption,maxheight=value
                with PILImage.open(OUT/filename) as im: w,h=im.size
                scale=min(width/w,maxheight/h)
                story.extend([Image(str(OUT/filename),width=w*scale,height=h*scale),Paragraph(escape(caption),styles['caption'])])
    def footer(canvas,doc):
        canvas.setStrokeColor(colors.HexColor('#b9cbd0')); canvas.line(40,35,A4[0]-40,35)
        canvas.setFont('CJK',8); canvas.setFillColor(colors.HexColor('#50626b'))
        canvas.drawString(40,23,'TIRT 光学验证 | 研究来源：闻建光团队 | 2026-10-24')
        canvas.drawRightString(A4[0]-40,23,str(doc.page))
    SimpleDocTemplate(str(PDF),pagesize=A4,leftMargin=40,rightMargin=40,topMargin=35,bottomMargin=47,
                      title='TIRT 光学验证报告',author='TIRT 项目').build(story,onFirstPage=footer,onLaterPages=footer)
    reader=PdfReader(PDF)
    assert len(reader.pages)==14,len(reader.pages)
    content='\n'.join(p.extract_text() for p in reader.pages)
    for token in ['闻建光团队','1,204','602','670','860','1e-6','0.084806']:
        assert token in content,token
    (OUT/'report_verification.json').write_text(json.dumps({
        'pdf_pages':len(reader.pages),'scene_charts':10,'figures':13,'optical_pairs':1204,
        'source_hashes_verified':True,'metrics_recalculated':True,
        'outputs':{p.relative_to(ROOT).as_posix():sha(p) for p in (MD,PDF)},
        'visual_review':'pending'
    },ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(f'Built {len(reader.pages)} pages; verified 1204 pairs and 20 scene-band metrics.')


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--charts',action='store_true')
    parser.add_argument('--pdf',action='store_true')
    args=parser.parse_args()
    scenes,paired,scores,totals,coverage=load()
    if args.charts or not args.pdf: charts(scenes,paired)
    if args.pdf or not args.charts: document(scenes,scores,totals,coverage)
