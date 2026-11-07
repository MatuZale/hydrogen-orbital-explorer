#!/usr/bin/env python3
"""Render C++ point clouds as a still and an optional camera-orbit movie.

Requires numpy and Pillow; movie output additionally uses imageio-ffmpeg.
Cloud geometry and wavefunction signs come only from export_cloud. Soft point
rendering and colours are presentation choices, not an isosurface calculation.
"""
from argparse import ArgumentParser
from pathlib import Path
import math
import subprocess

import numpy as np
from PIL import Image, ImageDraw, ImageFilter, ImageFont


def font(size, bold=False):
    candidates = [
        '/System/Library/Fonts/Supplemental/Arial Bold.ttf' if bold else '/System/Library/Fonts/Supplemental/Andale Mono.ttf',
        '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf' if bold else '/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf',
    ]
    for candidate in candidates:
        if Path(candidate).exists():
            return ImageFont.truetype(candidate, size)
    return ImageFont.load_default(size=size)


def cloud_layer(data, yaw, width, height):
    # Same fixed cloud each frame: only the viewing basis rotates.
    xyz = data[:, :3]
    cy, sy = math.cos(yaw), math.sin(yaw)
    cp, sp = math.cos(.35), math.sin(.35)
    right = np.array([cy, -sy, 0.0])
    forward = np.array([sy*cp, cy*cp, sp])
    up = np.cross(right, forward)
    scale = min(width, height) * .43 / np.quantile(np.linalg.norm(xyz, axis=1), .985)
    px = xyz @ right * scale + width*.5
    py = height*.5 - xyz @ up * scale
    x, y = np.floor(px).astype(int), np.floor(py).astype(int)
    ok = (x >= 0) & (x < width-1) & (y >= 0) & (y < height-1)
    x, y, px, py = x[ok], y[ok], px[ok], py[ok]
    positive = data[ok, 3] >= 0
    layers = np.zeros((2, height, width), np.float32)
    # Bilinear splatting removes pixel-grid popping during camera rotation.
    dx, dy = px-x, py-y
    for sx, sy, weight in [(0,0,(1-dx)*(1-dy)), (1,0,dx*(1-dy)), (0,1,(1-dx)*dy), (1,1,dx*dy)]:
        np.add.at(layers, (positive.astype(int), y+sy, x+sx), weight)
    negative, positive = layers
    # Saturating exposure; the visual is a projection of sampled probability.
    neg = 1-np.exp(-negative*.75)
    pos = 1-np.exp(-positive*.75)
    rgb = pos[:,:,None]*np.array([75.,177.,255.]) + neg[:,:,None]*np.array([255.,124.,55.])
    sharp = Image.fromarray(np.clip(rgb,0,255).astype(np.uint8))
    glow = np.asarray(sharp.filter(ImageFilter.GaussianBlur(2.0))).astype(float)
    background = np.array([6.,12.,24.])
    return Image.fromarray(np.clip(np.asarray(sharp)*.90+glow*.25+background,0,255).astype(np.uint8))


def frame(clouds, angle, size):
    w, h = size
    image = Image.new('RGB',size,(6,12,24))
    draw = ImageDraw.Draw(image)
    small, title, label = font(max(10,int(w*.012))), font(int(w*.029),True), font(int(w*.023),True)
    draw.text((int(w*.035),int(h*.045)), 'HYDROGEN / ORBITAL ATLAS',font=small,fill=(120,164,185))
    draw.text((int(w*.035),int(h*.09)), 'The shape of a wavefunction',font=title,fill=(237,241,230))
    top, bottom = int(h*.20),int(h*.79)
    panelw = w//3
    for i,(name,indices,data) in enumerate(clouds):
        image.paste(cloud_layer(data,angle+i*.12,panelw,bottom-top),(i*panelw,top))
        x = i*panelw+int(w*.035)
        draw.text((x,int(h*.80)),name,font=label,fill=(238,240,227))
        draw.text((x,int(h*.85)),indices,font=small,fill=(125,163,186))
    draw.line((int(w*.035),int(h*.91),int(w*.965),int(h*.91)), fill=(35,60,77),width=1)
    draw.text((int(w*.035),int(h*.94)),'BLUE + / AMBER -   |   REAL WAVEFUNCTION SIGN',font=small,fill=(150,184,199))
    return image


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--exporter',type=Path,default=Path('build/export_cloud'))
    parser.add_argument('--output',type=Path,default=Path('media'))
    parser.add_argument('--movie',action='store_true')
    parser.add_argument('--seconds',type=int,default=12)
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    cache=args.output/'clouds'
    cache.mkdir(exist_ok=True)
    clouds=[]
    for name,n,l,m in [('2p_x',2,1,1),('3d_z2',3,2,0),('4f_xyz',4,3,-2)]:
        path=cache/f'{name}.csv'
        subprocess.run([str(args.exporter.resolve()),str(path),str(n),str(l),str(m),'70000','42'],check=True)
        clouds.append((name,f'n={n} / l={l} / real basis',np.loadtxt(path,delimiter=',',skiprows=1)))
    frame(clouds,.60,(1280,720)).save(args.output/'orbitals-preview.webp',quality=92,method=6)
    if args.movie:
        import imageio_ffmpeg
        fps=30
        writer=imageio_ffmpeg.write_frames(str(args.output/'orbitals-rotation.mp4'),(1024,576),fps=fps,codec='libx264',quality=5,macro_block_size=16,output_params=['-movflags','+faststart'])
        writer.send(None)
        for i in range(args.seconds*fps):
            writer.send(np.asarray(frame(clouds,.60+2*math.pi*i/(args.seconds*fps),(1024,576))))
        writer.close()
    print(f'Rendered previews to {args.output}')


if __name__=='__main__':
    main()
