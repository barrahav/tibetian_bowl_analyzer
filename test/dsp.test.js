import test from 'node:test';
import assert from 'node:assert/strict';
import { analyze, preprocess } from '../src/dsp.js';
import { notePosition, noteLabel } from '../src/notes.js';

function signal(sr=16000,duration=5,parts=[[440,0.35,1.8],[660,0.18,1.2],[880,0.1,2.4]]){const n=Math.round(sr*duration),x=new Float32Array(n);for(let i=0;i<n;i++){const t=i/sr;for(const [f,a,tau] of parts)x[i]+=a*Math.exp(-t/tau)*Math.sin(2*Math.PI*f*t)}return x}
test('removes DC and reports clipping',()=>{const x=signal();for(let i=0;i<x.length;i++)x[i]+=0.1;const p=preprocess([x],16000);assert.ok(Math.abs(p.x.reduce((a,b)=>a+b,0)/x.length)<1e-8);assert.equal(p.clippedFraction,0)});
test('finds stable spectral components and independent decay estimates',()=>{const r=analyze([signal()],16000);assert.ok(r.components.length>=2);const freqs=r.components.map(c=>c.frequency);assert.ok(freqs.some(f=>Math.abs(f-440)<5));for(const c of r.components){assert.ok(c.decay.envelope.length>0);assert.ok(c.decay.fit);assert.ok(Number.isFinite(c.decay.fit.slopeDbPerSecond))}const slopes=r.components.map(c=>c.decay.fit.slopeDbPerSecond);assert.ok(new Set(slopes.map(t=>t.toFixed(1))).size>1,'each component should have an independent decay estimate')});
test('does not claim a natural decay for continuous rim excitation',()=>{const r=analyze([signal()],16000,'rim');assert.ok(r.components.length);assert.ok(r.components.every(c=>c.decay.fit===null));assert.match(r.components[0].decay.status,/Sustained excitation/)});
test('rejects quiet and short recordings',()=>{assert.throws(()=>analyze([new Float32Array(500)],16000),/too short/);assert.throws(()=>analyze([new Float32Array(16000)],16000),/too quiet/)});
test('places pitch between adjacent semitones on a logarithmic scale',()=>{const a=notePosition(440),middle=notePosition(Math.sqrt(440*466.1637615));assert.equal(noteLabel(a.nearest),'A4');assert.equal(a.nearestHz,440);assert.equal(noteLabel(a.lower),'A4');assert.equal(noteLabel(a.upper),'B♭4');assert.ok(Math.abs(a.position)<1e-7);assert.ok(Math.abs(middle.position-50)<0.001);assert.ok(Math.abs(notePosition(466.1637615).position-100)<0.001)});
