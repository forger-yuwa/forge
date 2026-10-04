# Salome スパイク: 隅ごとに独立な r(x) のフィレットを「ロフトしたガセット固体」で作り、流体 = 箱 − 固体 を
# NETGEN tet + ViscousLayers (prism) で切る。上隅 = 徐変で始まり終端まで付けたまま (鈍頭)。下隅 = 一定 → 終端で消す。
import sys, json, math
import salome; salome.salome_init()
from salome.geom import geomBuilder; from salome.smesh import smeshBuilder
import SMESH
gp=geomBuilder.New(); sm=smeshBuilder.New()
OUT=sys.argv[1] if len(sys.argv)>1 else "/tmp/junction"
H=0.1; Lsw,Lcowl,Lx=1.2*H,1.6*H,2.4*H; hw=1.0*H; tsw=0.05*H; tc=0.02*H
V=gp.MakeVertex; box=lambda a,b: gp.MakeBoxTwoPnt(V(*a),V(*b))
ramp=box((0,1.0*H,0),(Lx,1.3*H,2.0*H)); cowl=box((0,-tc,0),(Lcowl,0,hw+tsw)); sw=box((0,0,hw),(Lsw,1.0*H,hw+tsw))
def gusset(ys, sgn, prof):
    """隅 S=(ys, hw) のガセット固体。sgn=-1: 上隅 (流体は −y 側), +1: 下隅 (流体は +y 側)。prof = [(x, r), ...]"""
    secs=[]
    for x,r in prof:
        Sv=V(x,ys,hw); A=V(x,ys,hw-r); B=V(x,ys+sgn*r,hw)
        c=(ys+sgn*r, hw-r); M=V(x, c[0]-sgn*r*math.sin(math.pi/4), c[1]+r*math.cos(math.pi/4))   # 円弧中点 (S 側へ膨らむ)
        w=gp.MakeWire([gp.MakeLineTwoPnt(Sv,A), gp.MakeArc(A,M,B), gp.MakeLineTwoPnt(B,Sv)])
        secs.append(w)
    return gp.MakeThruSections(secs, True, 1e-7, False)
up=gusset(1.0*H,-1,[(0.0,0.01*H),(0.4*H,0.10*H),(0.8*H,0.10*H),(Lsw,0.10*H)])      # 徐変で始まり、終端まで付けたまま
lo=gusset(0.0,+1,[(0.0,0.08*H),(0.4*H,0.08*H),(0.8*H,0.08*H),(Lsw,0.01*H)])        # 一定 → 終端で (ほぼ) 消す
solid=gp.MakeFuseList([ramp,cowl,sw,up,lo],True,True)
fluid=gp.MakeCutList(box((0,-1.0*H,0),(Lx,1.3*H,2.0*H)),[solid],True)
props=gp.BasicProperties(fluid); nf=len(gp.SubShapeAll(fluid,gp.ShapeType["FACE"])); ns=len(gp.SubShapeAll(fluid,gp.ShapeType["SOLID"]))
print("GEOM: 流体 solid %d 個, 面 %d 枚, 体積 %.6e m3, valid=%s"%(max(ns,1),nf,props[2],gp.CheckShape(fluid)))
gp.ExportSTEP(fluid, OUT+".step")
# ---- メッシュ: NETGEN 1D-2D-3D + ViscousLayers (壁 = 箱の外面以外の面) ----
faces=gp.SubShapeAll(fluid,gp.ShapeType["FACE"])
def is_outer(f):
    (x0,x1,y0,y1,z0,z1)=gp.BoundingBox(f); e=1e-5
    return (abs(x0-x1)<e and (abs(x0)<e or abs(x0-Lx)<e)) or (abs(y0-y1)<e and abs(y0+1.0*H)<e) or (abs(z0-z1)<e and (abs(z0)<e or abs(z0-2.0*H)<e))
walls=[f for f in faces if not is_outer(f)]; outer=[f for f in faces if is_outer(f)]
print("壁面 %d 枚 / 外枠面 %d 枚"%(len(walls),len(outer)))
cfg=json.load(open(sys.argv[2])) if len(sys.argv)>2 else {}
h1=cfg.get("h1",2.0e-5); nl=cfg.get("nl",20); st=cfg.get("stretch",1.2)
tot=h1*(st**nl-1)/(st-1)
mesh=sm.Mesh(fluid); ng=mesh.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
p=ng.Parameters(); p.SetMaxSize(cfg.get("maxh",0.08*H)); p.SetMinSize(cfg.get("minh",0.004*H)); p.SetFineness(3); p.SetOptimize(1)
if cfg.get('wall_h'):
    gp.addToStudy(fluid,'fluid')
    for i,f in enumerate(walls):
        gp.addToStudyInFather(fluid,f,'wall_%d'%i); p.SetLocalSizeOnShape(f, cfg['wall_h'])
ids=[gp.GetSubShapeID(fluid,f) for f in outer]
vl=ng.ViscousLayers(tot, nl, st, ids, True)       # True = ids は「層を張らない面」
ok=mesh.Compute(); print("Compute:",ok)
info=mesh.GetMeshInfo(); cnt={str(k):v for k,v in info.items() if v}
print("要素:",{k.replace("Entity_",""):v for k,v in cnt.items() if any(s in k for s in ("Tetra","Penta","Pyramid","Hexa","Node"))})
print("VL 指定: 第一層 %.2e m × %d 層 (成長率 %.2f) = 総厚 %.3e m"%(h1,nl,st,tot))
try:
    errs=mesh.GetComputeErrors()
    for e in errs[:5]: print("  ERR:",e.comment[:150])
except Exception as ex: print("(errors n/a)",ex)
def fgroup(name, fs):
    if not fs: return
    gg=gp.CreateGroup(fluid, gp.ShapeType["FACE"]); gp.UnionList(gg, fs); gp.addToStudyInFather(fluid, gg, name)
    mesh.GroupOnGeom(gg, name, SMESH.FACE)
def cls(f):
    (x0,x1,y0,y1,z0,z1)=gp.BoundingBox(f); e=1e-5
    if abs(x0-x1)<e and abs(x0)<e: return "inlet"
    if abs(x0-x1)<e and abs(x0-Lx)<e: return "outlet"
    if abs(z0-z1)<e and abs(z0)<e: return "symmetry"
    return "plume_far"
if not cfg.get('wall_h'): gp.addToStudy(fluid,'fluid')
by={}
for f in outer: by.setdefault(cls(f),[]).append(f)
for k,v in by.items(): fgroup(k,v)
fgroup("wall", walls)
mesh.GroupOnGeom(fluid, "fluid", SMESH.VOLUME)
mesh.ExportMED(OUT+".med",auto_groups=False)
print("MED 書き出し:",OUT+".med")
