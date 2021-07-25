//dan olson
//2020

const frag_buffer = `

varying vec2 uVu;
uniform sampler2D tex;
uniform float time;

#define STEPS 255
#define EPS 0.00001
#define NEAR 0.
#define FAR 75.;

mat2 rot(float a) {

    float c = cos(a);
    float s = sin(a);
    
    return mat2(c,-s,s,c);
}

mat3 camOrthographic(vec3 ro,vec3 ta,float r) {
     
     vec3 w = normalize(ta - ro); 
     vec3 p = vec3(sin(r),cos(r),0.);           
     vec3 u = normalize(cross(w,p)); 
     vec3 v = normalize(cross(u,w));

     return mat3(u,v,w); 
} 

vec2 opu(vec2 d1,vec2 d2) {
    return (d1.x < d2.x) ? d1 : d2;
} 

float box(vec3 p,vec3 b) {

    vec3 d = abs(p) - b;
    return length(max(d,0.0)) + min(max(d.x,max(d.y,d.z)),0.0);
}


vec2 scene(vec3 p) {

    vec2 res = vec2(1.,0.);

    float d = 0.;     
    d = box(p,vec3(1.));

    res = opu(res,vec2(d,2.)); 
    res = opu(res,vec2(-box(p,vec3(25.)),3.));

  return res;

}

vec2 rayScene(vec3 ro,vec3 rd) {
    
    float d = -1.0;
    float s = NEAR;
    float e = FAR;  

    for(int i = 0; i < STEPS; i++) {

        vec3 p = ro + s * rd;
        vec2 dist = scene(p);
   
        if(abs(dist.x) < EPS || e <  dist.x ) { break; }
        s += dist.x;
        d = dist.y;

        }
 
        if(e < s) { d = -1.0; }
        return vec2(s,d);

}

vec3 calcNormal(vec3 p) {

    vec2 e = vec2(1.,-1.)*EPS;

    return normalize(vec3(
    vec3(e.x,e.y,e.y) * scene(p + vec3(e.x,e.y,e.y)).x +
    vec3(e.y,e.x,e.y) * scene(p + vec3(e.y,e.x,e.y)).x +
    vec3(e.y,e.y,e.x) * scene(p + vec3(e.y,e.y,e.x)).x + 
    vec3(e.x,e.x,e.x) * scene(p + vec3(e.x,e.x,e.x)).x

    ));
    
}

vec3 renderScene(vec3 ro,vec3 rd) {
 
vec2 d = rayScene(ro, rd);

vec3 col = vec3(1.) * max(0.,rd.y);

if(d.y >= 0.) {

vec3 p = ro + rd * d.x;
vec3 n = calcNormal(p);
vec3 l = normalize(vec3(10.));

float dif = dot(n,l)*0.5+0.5;

vec3 col_xy = texture(tex,p.xy*0.5+0.5).rgb;
vec3 col_xz = texture(tex,p.xz*0.5+0.5).rgb;
vec3 col_zy = texture(tex,p.zy*0.5+0.5).rgb;

n = abs(n);

col = col_xy*n.z + col_xz*n.y + col_zy*n.x;

col *= dif;
col = mix(col,vec3(1.),1.-exp(-0.00001 * d.x*d.x*d.x)); 

}

return col;
}

void main() {

vec3 color = vec3(0.);

vec3 ro = vec3(2.);
ro.xy *= rot(time*0.15);

vec3 ta = vec3(0.0);

vec2 uv = uVu.xy * 2. - 1.; 

mat3 cm = camOrthographic(ro,ta,0.);
vec3 rd = cm * normalize(vec3(uv,2.));

vec3 col = renderScene(ro,rd);

col = pow(col,vec3(0.4545));
color += col;

gl_FragColor = vec4(color,1.0);
}

`;

const frag = `

uniform vec2 resolution;
uniform sampler2D tex;
uniform float time;

#define STEPS 255
#define EPS 0.00001
#define NEAR 0.
#define FAR 500.
#define AA 2
#define seed 5425122

float h11(float p) {
    uvec2 n = uint(int(p)) * uvec2(uint(int(seed)),2531151992.0);
    uint h = (n.x ^ n.y) * uint(int(seed));
    return float(h) * (1./float(0xffffffffU));
}

float h21(vec2 p) {
    uvec2 n = uvec2(ivec2(p)) * uvec2(uint(int(seed)),2531151992.0);
    uint h = (n.x ^ n.y) * uint(int(seed));
    return float(h) * (1./float(0xffffffffU));
}

float n3(vec3 x) {
    vec3 p = floor(x);
    vec3 f = fract(x);

    f = f * f * (3.0 - 2.0 * f);
    float n = p.x + p.y * 157.0 + 113.0 * p.z;

    return mix(mix(mix(h11(n + 0.0), 
                       h11(n + 1.0),f.x),
                   mix(h11(n + 157.0),
                       h11(n + 158.0),f.x),f.y),
               mix(mix(h11(n + 113.0), 
                       h11(n + 114.0),f.x),
                   mix(h11(n + 270.0), 
                       h11(n + 271.0),f.x),f.y),f.z);
}

float f3(vec3 x,int octaves,float hurst) {
    float s = 0.;
    float h = exp2(-hurst);
    float f = 1.;
    float a = .5;

    for(int i = 0; i < octaves; i++) {

        s += a * n3(f * x);  
        f *= 2.;
        a *= h;
    }
    return s;
}

mat2 rot(float a) {

    float c = cos(a);
    float s = sin(a);
    
    return mat2(c,-s,s,c);
}

mat3 camOrthographic(vec3 ro,vec3 ta,float r) {
     
     vec3 w = normalize(ta - ro); 
     vec3 p = vec3(sin(r),cos(r),0.);           
     vec3 u = normalize(cross(w,p)); 
     vec3 v = normalize(cross(u,w));

     return mat3(u,v,w); 
} 

vec2 opu(vec2 d1,vec2 d2) {
    return (d1.x < d2.x) ? d1 : d2;
} 

float sphere(vec3 p,float r) { 
     
    return length(p) - r;
}

float box(vec3 p,vec3 b) {

    vec3 d = abs(p) - b;
    return length(max(d,0.0)) + min(max(d.x,max(d.y,d.z)),0.0);
}

vec2 scene(vec3 p) {

    vec2 res = vec2(1.,0.);

    float d = 0.;     
    d = -box(p,vec3(10.));
    res = opu(res,vec2(d,10.)); 
    res = opu(res,vec2(box(p,vec3(1.)),2.)); 
    

    return res;

}

vec2 rayScene(vec3 ro,vec3 rd) {
    
    float d = -1.0;
    float s = NEAR;
    float e = FAR;  

    for(int i = 0; i < STEPS; i++) {

        vec3 p = ro + s * rd;
        vec2 dist = scene(p);
   
        if(abs(dist.x) < EPS || e <  dist.x ) { break; }
        s += dist.x;
        d = dist.y;

        }
 
        if(e < s) { d = -1.0; }
        return vec2(s,d);

}

float shadow(vec3 ro,vec3 rd) {

    float res = 1.0;
    float t = 0.005;
    float ph = 1e10;
    
    for(int i = 0; i < 45; i++ ) {
        
        float h = scene(ro + rd * t  ).x;

        float y = h * h / (2. * ph);
        float d = sqrt(h*h-y*y);         
        res = min(res,25. * d/max(0.,t-y));
        ph = h;
        t += h;
    
        if(res < EPS || t > 10.) { break; }

        }


        return clamp(res,0.0,1.0);

}

vec3 calcNormal(vec3 p) {

    vec2 e = vec2(1.,-1.)*EPS;

    return normalize(vec3(
    vec3(e.x,e.y,e.y) * scene(p + vec3(e.x,e.y,e.y)).x +
    vec3(e.y,e.x,e.y) * scene(p + vec3(e.y,e.x,e.y)).x +
    vec3(e.y,e.y,e.x) * scene(p + vec3(e.y,e.y,e.x)).x + 
    vec3(e.x,e.x,e.x) * scene(p + vec3(e.x,e.x,e.x)).x

    ));
    
}

vec3 renderScene(vec3 ro,vec3 rd) {
 
vec2 d = rayScene(ro, rd);

vec3 col = vec3(1.) * max(0.,rd.y);

if(d.y >= 0.) {

vec3 p = ro + rd * d.x;
vec3 n = calcNormal(p);

vec3 l = normalize(vec3(3.));

vec3 h = normalize(l - rd);
vec3 r = reflect(rd,n);

float amb = clamp(0.5 + 0.5 * n.y,0.,1.);

float dif = clamp(dot(n,l),0.0,1.0);

float spe = pow(clamp(dot(n,h),0.0,1.0),16.)
* dif * (.04 + 0.9 * pow(clamp(1. + dot(h,rd),0.,1.),5.));

float fre = pow(clamp(1. + dot(n,rd),0.0,1.0),2.0);
float ref = smoothstep(-.2,.2,r.y);

vec3 linear = vec3(0.);

dif *= shadow(p,l);
ref *= shadow(p,r);

linear += dif * vec3(.5);
linear += amb * vec3(0.01,0.05,0.05);

if(d.y == 10.) {
    col = vec3(.5);
}

if(d.y == 2.) {

    float nl; 
    nl += f3(p+f3(p,4,h11(390.)),6,h11(35.));
    col += vec3(nl);

}

col = col * linear;
col += spe * vec3(1.,0.97,1.); 

col = mix(col,vec3(1.),1.-exp(-0.00001 * d.x*d.x*d.x)); 

}

return col;
}

void main() {

vec3 color = vec3(0.);

vec3 ro = vec3(3.);
vec3 ta = vec3(0.0);

ro.xz *= rot(time*0.1);

for(int k = 0; k < AA; ++k) {
    for(int l = 0; l < AA; ++l) {

        vec2 o = vec2(float(l),float(k)) / float(AA) - .5;

        vec2 uv = (2. * (gl_FragCoord.xy + o) -
        resolution.xy) / resolution.y; 

        mat3 cm = camOrthographic(ro,ta,0.);
        vec3 rd = cm * normalize(vec3(uv.xy,2.));
        vec3 col = renderScene(ro,rd);    

        col = pow(col,vec3(0.4545));
        color += col;
    }
}

color /= float(AA*AA);
gl_FragColor = vec4(color,1.0);
}

`;

const vert = `
 
varying vec2 uVu;

void main() {
 
uVu = uv;
gl_Position = vec4(position,1.);

}

`;

import {

WebGLRenderer,
Vector2,
Scene, 
PerspectiveCamera,
Mesh,
Clock,
PlaneBufferGeometry,
ShaderMaterial,
Texture,
WebGLRenderTarget

} from
'https://cdnjs.cloudflare.com/ajax/libs/three.js/r123/three.module.js';

let w,h; 

let scene;
let scene_rtt;

let cam;

let renderer;

let uniforms;

let material;
let material_texture;

let buffer;

let res;

let clock;

class Main { 

    constructor(container) {

        scene = new Scene();       
        scene_rtt = new Scene();
        
        cam = new PerspectiveCamera();
        renderer = new WebGLRenderer(); 
        renderer.setPixelRatio(window.devicePixelRatio);

        container.append(renderer.domElement);

        this.resize();

        clock = new Clock();

        const plane = new PlaneBufferGeometry(2,2);

        buffer = new WebGLRenderTarget(w,h,{
        minFilter : Texture.LinearFilter,
        magFilter : Texture.NearestFilter,
        format    : Texture.RGBFormat 
        });

        material = new ShaderMaterial({

            uniforms : {

                resolution : { 
                value : new Vector2(w,h)
                },

                time : { value : 0. }
            },

            vertexShader : vert,
            fragmentShader : frag

        });

        material.glslVersion = 'THREE.GLSL3';
 
        material_texture = new ShaderMaterial({

            uniforms : {

                tex  : { value : buffer.texture },
                time : { value : 0. }
            
            },
        
        vertexShader : vert,
        fragmentShader : frag_buffer

        });

        material_texture.glslVersion = 'THREE.GLSL3';

        const mesh_texture = new Mesh(plane,material_texture);
        scene.add(mesh_texture);

        const mesh = new Mesh(plane,material);
        scene_rtt.add(mesh);



    }

    begin() {

        renderer.setAnimationLoop(() => {

            renderer.setRenderTarget(buffer);
            renderer.clear();
            renderer.render(scene_rtt,cam);
            
            renderer.setRenderTarget(null);
            renderer.clear();
            renderer.render(scene,cam);

            this.resize();        
    
            material.uniforms.time.value = clock.getElapsedTime(); 

            material_texture.uniforms.time.value = clock.getElapsedTime();

            });

    }

    end() { 
        renderer.setAnimationLoop(null);
    }   

    render() {
        this.begin();
    }

    resize() {

        w = window.innerWidth;
        h = window.innerHeight;

        cam.aspect = w/h;
        renderer.setSize(w,h);     

    }

}

export { Main };

