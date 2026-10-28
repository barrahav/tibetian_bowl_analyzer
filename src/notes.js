const noteNames=['C','C♯','D','E♭','E','F','F♯','G','A♭','A','B♭','B'];

export function notePosition(frequency){
  if(!Number.isFinite(frequency)||frequency<=0)throw new Error('תדר חייב להיות מספר חיובי.');
  const midi=69+12*Math.log2(frequency/440),nearest=Math.round(midi),nearestHz=440*2**((nearest-69)/12),delta=1200*Math.log2(frequency/nearestHz);
  const lower=Math.floor(midi);
  const upper=lower+1,loHz=440*2**((lower-69)/12),hiHz=440*2**((upper-69)/12),position=100*Math.log(frequency/loHz)/Math.log(hiHz/loHz);
  return{lower,upper,loHz,hiHz,position:Math.max(0,Math.min(100,position)),nearest,nearestHz,delta};
}

export function noteLabel(midi){const pc=((midi%12)+12)%12,octave=Math.floor(midi/12)-1;return`${noteNames[pc]}${octave}`}
