// Local acoustic analysis. No mechanical-mode inference is made.
export const SNR_THRESHOLD_DB = 10;
export function preprocess(channels, sampleRate) {
  if (!channels?.length || channels[0].length < sampleRate * 0.5) throw new Error('Recording is too short (need at least 0.5 seconds).');
  const n = Math.min(...channels.map(c => c.length));
  const x = new Float64Array(n);
  for (const c of channels) for (let i=0;i<n;i++) x[i] += c[i] / channels.length;
  const mean = x.reduce((a,b)=>a+b,0)/n;
  let peak=0, clipped=0, sum=0;
  for(let i=0;i<n;i++){ x[i]-=mean; const a=Math.abs(x[i]); peak=Math.max(peak,a); sum+=x[i]*x[i]; if(a>=0.999)clipped++; }
  const rms=Math.sqrt(sum/n);
  if(rms<1e-5) throw new Error('Recording is too quiet to analyze reliably.');
  return {x,sampleRate,duration:n/sampleRate,peak,rms,clippedFraction:clipped/n, warnings:[...(clipped/n>0.001?['Possible clipping detected.']:[]),...(rms<0.001?['Low input level; weak components may be missed.']:[])]};
}
export function fft(re, im) {
  const n=re.length;
  if(n!==im.length || (n&(n-1))) throw new Error('FFT length must be a power of two.');
  for(let i=1,j=0;i<n;i++){let bit=n>>1;for(;j&bit;bit>>=1)j^=bit;j^=bit;if(i<j){[re[i],re[j]]=[re[j],re[i]];[im[i],im[j]]=[im[j],im[i]];}}
  for(let len=2;len<=n;len*=2){const ang=-2*Math.PI/len,wr0=Math.cos(ang),wi0=Math.sin(ang);for(let i=0;i<n;i+=len){let wr=1,wi=0;for(let j=0;j<len/2;j++){const u=i+j,v=u+len/2,tr=wr*re[v]-wi*im[v],ti=wr*im[v]+wi*re[v];re[v]=re[u]-tr;im[v]=im[u]-ti;re[u]+=tr;im[u]+=ti;const z=wr*wr0-wi*wi0;wi=wr*wi0+wi*wr0;wr=z;}}}
  return {re,im};
}
function spectrum(x, sr, start, size) {
  const re=new Float64Array(size), im=new Float64Array(size);
  for(let i=0;i<size;i++){const k=start+i, w=0.5-0.5*Math.cos(2*Math.PI*i/(size-1));re[i]=(x[k]||0)*w;}
  fft(re,im); return Array.from({length:size/2},(_,i)=>Math.hypot(re[i],im[i]));
}
function peaksFor(x,sr,start,size){
  const m=spectrum(x,sr,start,size), df=sr/size, peaks=[];
  for(let i=3;i<m.length-3;i++) if(m[i]>m[i-1]&&m[i]>=m[i+1]){
    const local=Math.max(1e-12, [...m.slice(i-12,i-3),...m.slice(i+4,i+13)].sort((a,b)=>a-b)[8]||1e-12);
    const snrDb=20*Math.log10(Math.max(m[i],1e-12)/local);
    if(snrDb<SNR_THRESHOLD_DB)continue;
    const d=0.5*(Math.log(Math.max(m[i-1],1e-12))-Math.log(Math.max(m[i+1],1e-12)))/(Math.log(Math.max(m[i-1],1e-12))-2*Math.log(Math.max(m[i],1e-12))+Math.log(Math.max(m[i+1],1e-12))||1);
    peaks.push({frequency:(i+d)*df, amplitude:m[i], prominence:snrDb, snrDb, binSpacing:df});
  }
  return peaks.sort((a,b)=>b.amplitude-a.amplitude);
}
function envelope(x,sr,f,hop){
  const radius=Math.max(128,Math.round(sr*0.06)), times=[], amps=[];
  // Narrowband complex demodulation with a Hann low-pass window; independent envelope per component.
  for(let c=radius;c<x.length-radius;c+=hop){let re=0,im=0,wSum=0;for(let j=-radius;j<=radius;j++){const w=0.5+0.5*Math.cos(Math.PI*j/radius),v=x[c+j]*w,ph=2*Math.PI*f*j/sr;re+=v*Math.cos(ph);im-=v*Math.sin(ph);wSum+=w;}times.push(c/sr);amps.push(2*Math.hypot(re,im)/wSum);}
  return {times,amps};
}
function fitDecay(env,noiseAmp){
  const {times,amps}=env, max=Math.max(...amps), refIdx=amps.findIndex(a=>a>=max*0.8), ref=Math.max(amps[refIdx]||max,1e-12), floor=Math.max(noiseAmp,ref*1e-5), levels=amps.map(a=>20*Math.log10(Math.max(a,floor)/ref));
  let end=levels.findIndex((v,i)=>i>refIdx&&v<=-25); if(end<0)end=levels.findIndex((v,i)=>i>refIdx&&amps[i]<=noiseAmp*2); if(end<0)end=levels.length-1;
  const begin=Math.min(levels.length-1,refIdx+2); let lo=begin,hi=end;
  // Fit a common early decay region, truncated at the estimated noise floor.
  if(hi-lo<8)return {reference:ref,envelope:amps,relativeDb:levels,fit:null,status:'Insufficient data for decay fit'};
  let sx=0,sy=0,sxx=0,sxy=0,n=0;for(let i=lo;i<=hi;i++){const t=times[i]-times[lo],y=levels[i];sx+=t;sy+=y;sxx+=t*t;sxy+=t*y;n++;}
  const slope=(n*sxy-sx*sy)/(n*sxx-sx*sx), intercept=(sy-slope*sx)/n;let ss=0,tot=0,mean=sy/n;for(let i=lo;i<=hi;i++){const t=times[i]-times[lo],y=levels[i];ss+=(y-(intercept+slope*t))**2;tot+=(y-mean)**2;}
  const r2=tot?1-ss/tot:0, range=levels[lo]-levels[hi];
  const fit={start:times[lo],end:times[hi],slopeDbPerSecond:slope,r2,rangeDb:range,t20:slope<0&&range>=18? -20/slope:null,t30:slope<0&&range>=27? -30/slope:null,t60:slope<0&&range>=27? -60/slope:null,points:times.slice(lo,hi+1).map((t,i)=>({t,db:intercept+slope*(t-times[lo])}))};
  const status=r2>=0.85&&range>=12?'Estimate with limitations':'Decay fit is uncertain'; return {reference:ref,envelope:amps,relativeDb:levels,fit,status};
}
export function analyze(channels,sampleRate,mode='soft'){
  const p=preprocess(channels,sampleRate), x=p.x, size=Math.min(16384,2**Math.floor(Math.log2(x.length))), start=Math.max(0,Math.floor((x.length-size)*0.35));
  const raw=peaksFor(x,sampleRate,start,size), nyquist=sampleRate/2;
  const candidates=[];
  for(const pk of raw){if(pk.frequency<35||pk.frequency>Math.min(8000,nyquist)||candidates.some(q=>Math.abs(q.frequency-pk.frequency)<Math.max(12,pk.frequency*.012)))continue;
    let persistence=0, deviations=[]; for(let w=0;w<4;w++){const at=Math.floor(w*(x.length-size)/3), local=peaksFor(x,sampleRate,at,size).filter(z=>Math.abs(z.frequency-pk.frequency)<Math.max(15,pk.frequency*.015)); if(local.length){persistence++;deviations.push(local[0].frequency);}}
    if(persistence<2)continue;
    const freq=deviations.reduce((a,b)=>a+b,0)/deviations.length;
    candidates.push({...pk,frequency:freq,persistence,score:pk.prominence+5*persistence});
  }
  candidates.sort((a,b)=>b.score-a.score);
  const selected=candidates.slice(0,3).sort((a,b)=>a.frequency-b.frequency), hop=Math.max(128,Math.round(sampleRate*0.02));
  const noiseAmp=Math.max(1e-7,p.rms*0.0002);
  const components=selected.map((pk,i)=>{const env=envelope(x,sampleRate,pk.frequency,hop),estimated=fitDecay(env,noiseAmp),decay=mode==='rim'?{...estimated,fit:null,status:'Sustained excitation; no natural decay inferred'}:estimated;return {id:String.fromCharCode(65+i),...pk,decay};});
  const dominant=components.reduce((a,b)=>!a||b.amplitude>a.amplitude?b:a,null), harmonic=components.length>=3&&components.slice(1).every((p,i)=>Math.abs(p.frequency/components[0].frequency-Math.round(p.frequency/components[0].frequency))<0.035);
  const bins=spectrum(x,sampleRate,start,size), spec=bins.map((a,i)=>({f:i*sampleRate/size,db:20*Math.log10(Math.max(a,1e-9))})).filter((_,i)=>i%2===0);
  return {...p,mode,fftSize:size,binSpacing:sampleRate/size,snrThresholdDb:SNR_THRESHOLD_DB,components,dominant:dominant?.id||null,lowest:components[0]?.id||null,candidateFundamental:harmonic?components[0].id:null, fundamentalStatus:harmonic?'Candidate; harmonic evidence present':'Undetermined; no sufficient harmonic evidence',spectrum:spec,quality:components.length?'Estimate with limitations':'Measurement failed: no stable tonal component exceeded the SNR threshold', warnings:[...p.warnings,...(mode==='rim'?['Rim playing remains sustained excitation: no natural decay is fitted during rubbing. Record the release clearly for any post-rub decay estimate.']:[])]};
}
