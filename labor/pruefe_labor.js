/* pruefe_labor.js -- Knopfprobe im echten Browser.
 *
 * Aufruf:  NODE_PATH=/tmp/node_modules node pruefe_labor.js [datei]
 *
 * Geprueft wird nicht, ob die Seite laedt, sondern ob sie TUT: jeder Reiter wird
 * geoeffnet, jeder Knopf gedrueckt, jeder Schieber bewegt, und die Ergebnisse werden
 * gegen harte Sollwerte gehalten. „Laedt ohne Fehler" ist keine Pruefung.
 */
const { chromium } = require('playwright');
const path = require('path');
const fs = require('fs');

/* Befund 01.11.2026: Hier stand fest „CFD_Labor_1.0.html". Wer die Probe ohne Argument
 * startete, pruefte eine ZWEI Fassungen alte Datei und bekam ein gruenes Ergebnis fuer
 * etwas, das gar nicht ausgeliefert wird. Eine Fassungsnummer gehoert an EINE Stelle --
 * hier wird stattdessen die hoechste vorhandene genommen. */
function neuesteFassung() {
  const w = path.join(__dirname, '..');
  const d = fs.readdirSync(w).filter(f => /^CFD_Labor_\d+\.\d+\.html$/.test(f));
  if (!d.length) throw new Error('keine CFD_Labor_*.html gefunden');
  d.sort((a, b) => {
    const z = s => s.match(/(\d+)\.(\d+)/).slice(1, 3).map(Number);
    const [a1, a2] = z(a), [b1, b2] = z(b);
    return a1 - b1 || a2 - b2;
  });
  return path.join(w, d[d.length - 1]);
}
const DATEI = process.argv[2] || neuesteFassung();

/* Die SOLLLISTE der Reiter -- aus den Aussagen des Auftraggebers, nicht aus der Seite.
 * Befund B17: Die Probe ging bisher `window.PRUEFSTAND.reiter` durch, also eine Liste
 * AUS DEM PRUEFLING. Damit kann sie pruefen, ob die vorhandenen Reiter funktionieren,
 * nie aber, ob ein geforderter fehlt -- die Kurzfassung fehlte in 1.1 bei 36/36 gruen.
 * Eine Soll-Liste muss ausserhalb des Prueflings stehen. */
const SOLL_REITER = ['labor', 'zeitverlauf', 'wirbel', 'kurz', 'verfahren', 'befunde', 'anleitung', 'quellen', 'lizenz'];
const BILD = process.env.BILDORDNER || '/tmp';

(async () => {
  const urteile = [];
  const U = (text, ok, beleg) => urteile.push([text, !!ok, beleg]);
  const browser = await chromium.launch();
  const seite = await browser.newPage({ viewport: { width: 1500, height: 1000 } });
  const fehler = [];
  seite.on('console', m => { if (m.type() === 'error') fehler.push(m.text()); });
  seite.on('pageerror', e => fehler.push('pageerror: ' + e.message));

  await seite.goto('file://' + path.resolve(DATEI));
  await seite.waitForTimeout(900);
  U('Seite laedt ohne Fehler in der Konsole', fehler.length === 0, fehler.slice(0, 2).join(' | ') || 'keine');

  const ps = await seite.evaluate(() => window.PRUEFSTAND);
  U('Diagnoseschnittstelle PRUEFSTAND vorhanden', !!ps && !!ps.fassung, ps ? ps.fassung : '—');

  // ---- die Reiterleiste gegen die SOLLLISTE halten (Befund B17) --------------------
  {
    const ist = (ps && ps.reiter) || [];
    const fehlt = SOLL_REITER.filter(r => !ist.includes(r));
    const zuviel = ist.filter(r => !SOLL_REITER.includes(r));
    U('alle geforderten Reiter sind da (Soll gegen Ist)', !fehlt.length && !zuviel.length,
      fehlt.length ? 'fehlt: ' + fehlt.join(', ')
                   : (zuviel.length ? 'unerwartet: ' + zuviel.join(', ')
                                    : ist.length + ' Reiter, wie gefordert'));
  }

  // ---- jeder Reiter laesst sich oeffnen und zeigt Inhalt ---------------------------
  for (const r of ps.reiter) {
    await seite.click('#r-' + r);
    await seite.waitForTimeout(260);
    const z = await seite.evaluate(n => {
      const e = document.getElementById('ansicht-' + n);
      return { da: !e.hidden, text: (e.innerText || '').trim().length };
    }, r);
    U(`Reiter „${r}" oeffnet und zeigt Inhalt`, z.da && z.text > 200, `${z.text} Zeichen`);
  }

  // ---- Kopf: Verfasser und Entstehung, und der Lizenzreiter (K17) ------------------
  {
    const kopf = await seite.evaluate(() => {
      const m = document.querySelector('.marke');
      const v = document.querySelector('.marke .verfasser');
      const r = v ? v.getBoundingClientRect() : null;
      return { text: m ? m.innerText.replace(/\s+/g, ' ') : '',
               sichtbar: !!r && r.height > 5 && r.width > 50 };
    });
    U('der Kopf nennt den Verfasser mit Titel',
      /Prof\. Dr\.-Ing\. Ralph Wystup/.test(kopf.text) && kopf.sichtbar,
      kopf.text.slice(0, 110));
    U('der Kopf nennt, dass es mit KI und Agent entstanden ist',
      /KI und Agent/.test(kopf.text) && /Claude Code/.test(kopf.text),
      kopf.text.includes('Claude Code') ? 'Claude Code, Anthropic genannt' : 'fehlt');

    await seite.click('#r-lizenz'); await seite.waitForTimeout(250);
    const liz = await seite.evaluate(() => {
      const e = document.getElementById('ansicht-lizenz');
      return { da: !e.hidden, text: (e.innerText || '').replace(/\s+/g, ' ') };
    });
    U('der Reiter „Lizenz" öffnet und trägt die Bedingungen des Urhebers',
      liz.da && /may not be sold/.test(liz.text) &&
      /include an acknowledgment/.test(liz.text) && /Shewchuk/.test(liz.text),
      `${liz.text.length} Zeichen`);
    U('der Lizenzreiter sagt, dass die eigene Freigabe keine Triangle-Lizenz erteilt',
      liz.da && /nicht.{0,30}erteil/i.test(liz.text),
      liz.da ? 'Satz vorhanden' : '—');
    U('der Lizenzreiter nennt beide Arbeiten zum Zitieren',
      liz.da && /Engineering a 2D Quality Mesh Generator/.test(liz.text) &&
      /Delaunay Refinement Algorithms/.test(liz.text),
      'beide Zitate');
    await seite.click('#r-labor'); await seite.waitForTimeout(200);
    await seite.click('#wahl-vergleich button:nth-child(5)'); await seite.waitForTimeout(250);
  }

  // ---- Programme: liegen die Quellen wirklich bei? (K18, Einwand Fred) ------------
  {
    await seite.click('#r-quellen'); await seite.waitForTimeout(300);
    const q = await seite.evaluate(() => {
      const s = document.getElementById('s-quelle');
      const namen = s ? [...s.options].map(o => o.textContent.split(' — ')[0]) : [];
      return { anzahl: namen.length, namen: namen,
               text: (document.getElementById('quelltext').textContent || '').length,
               kopf: document.getElementById('h-quelle-kopf').textContent };
    });
    /* Die Liste muss die Programme enthalten, auf die das Prüfprotokoll verweist.
     * „Irgendwelche Dateien sind da" genügt nicht -- die GENANNTEN müssen da sein. */
    const PFLICHT = ['stroemung/sg_stempel.py', 'stroemung/ns_loeser.py',
                     'stroemung/impuls_matrix.py', 'stroemung/wirbelstrasse.py',
                     'stroemung/netz_fem.py', 'stroemung/fem_galerkin.py',
                     'stroemung/fem_wirbelstrasse.py', 'stroemung/pruefprotokoll.py',
                     'labor/kern_cfd.js', 'labor/kern_fem.js', 'labor/netz_browser.js',
                     'labor/pruefe_kern_fem.js'];
    const fehlt = PFLICHT.filter(n => !q.namen.includes(n));
    U('die Programme liegen im Wortlaut bei, nicht als Verweis',
      q.anzahl >= 30 && !fehlt.length,
      fehlt.length ? 'fehlt: ' + fehlt.join(', ') : `${q.anzahl} Dateien`);
    U('der angezeigte Quelltext ist wirklich Quelltext', q.text > 1000,
      `${q.text} Zeichen`);
    /* Befund 01.11.2026: Die Programme lagen zwar alle bei, aber die Seite zeigte
     * immer nur EINES -- in einer Auswahlliste, die man aufklappen musste. Der Leser
     * musste denken, es liege nur ein Programm bei, und hat genau das gefragt.
     * „Vollständig" muss man SEHEN, nicht aufklappen. Deshalb: Die Übersicht muss
     * jede einzelne Datei als eigene Zeile führen, mit ihrer Aufgabe. */
    const tafel = await seite.evaluate(() => {
      const z = [...document.querySelectorAll('#quellliste tr[data-i]')];
      return {
        zeilen: z.length,
        ohneAufgabe: z.filter(r => (r.cells[1].textContent || '').trim().length < 10).length,
        ohneSumme: z.filter(r => !/^[0-9a-f]{12}/.test(r.cells[4].textContent || '')).length,
        kopf: document.getElementById('b-quellzahl').textContent,
      };
    });
    U('die Übersicht führt JEDE Datei als eigene Zeile, nicht nur die ausgewählte',
      tafel.zeilen === q.anzahl && tafel.zeilen >= 34,
      `${tafel.zeilen} Zeilen für ${q.anzahl} Dateien · ${tafel.kopf}`);
    U('zu jeder Datei steht ihre Aufgabe', tafel.ohneAufgabe === 0,
      tafel.ohneAufgabe ? `${tafel.ohneAufgabe} ohne Aufgabe` : 'alle beschrieben');
    U('zu jeder Datei steht ihre Prüfsumme in der Übersicht', tafel.ohneSumme === 0,
      tafel.ohneSumme ? `${tafel.ohneSumme} ohne Prüfsumme` : 'alle mit sha256');
    U('zu jeder Datei stehen Größe und Prüfsumme', /sha256 [0-9a-f]{64}/.test(q.kopf),
      q.kopf.split('\n')[0] || '—');
    // umschalten muss einen ANDEREN Text zeigen
    const zwei = await seite.evaluate(() => {
      const s = document.getElementById('s-quelle');
      const a = document.getElementById('quelltext').textContent;
      s.value = String(Math.min(5, s.options.length - 1));
      zeigeQuelle();
      return { a: a.length, b: document.getElementById('quelltext').textContent.length };
    });
    U('die Dateiauswahl schaltet wirklich um', zwei.a !== zwei.b,
      `${zwei.a} gegen ${zwei.b} Zeichen`);
    await seite.click('#r-anleitung'); await seite.waitForTimeout(250);
    const anl = await seite.evaluate(() =>
      (document.getElementById('ansicht-anleitung').innerText || '').replace(/\s+/g, ' '));
    U('das Blatt für Fred steht im Reiter, nicht nur im Verzeichnis',
      /Kurz für Fred/.test(anl) && /Beweisstück/.test(anl),
      `${anl.length} Zeichen im Anleitungsreiter`);
    await seite.click('#r-labor'); await seite.waitForTimeout(200);
    await seite.click('#wahl-vergleich button:nth-child(5)'); await seite.waitForTimeout(250);
  }

  // ---- Kurzfassung: die vier bestellten Abschnitte sind wirklich drin (K13) --------
  {
    await seite.click('#r-kurz');
    await seite.waitForTimeout(220);
    const k = await seite.evaluate(() => {
      const e = document.getElementById('ansicht-kurz');
      const ue = [...e.querySelectorAll('h1,h2,h3')].map(h => h.innerText.trim());
      return { ue, zeichen: (e.innerText || '').trim().length };
    });
    for (const soll of ['Aufgabe', 'Zielsetzung', 'Methode', 'Fazit']) {
      U(`Kurzfassung hat den Abschnitt „${soll}"`,
        k.ue.some(h => h.toLowerCase().startsWith(soll.toLowerCase())),
        k.ue.join(' · ').slice(0, 90));
    }
    U('Kurzfassung ist zwei Seiten lang, nicht ein Satz', k.zeichen > 3000,
      k.zeichen + ' Zeichen');
  }

  // ---- Vergleich 1: Zweigstempel ---------------------------------------------------
  await seite.click('#r-labor'); await seite.waitForTimeout(200);
  await seite.click('#wahl-vergleich button:nth-child(1)');
  await seite.waitForTimeout(300);
  for (const [sp, name] of [[20, 'Pe = 100'], [43, 'Pe ~ 20000'], [60, 'Pe = 1e6'], [10, 'Pe = 10']]) {
    await seite.evaluate(v => { document.getElementById('s-pe').value = String(v); rechneStempel(); }, sp);
    await seite.waitForTimeout(140);
    const t = await seite.evaluate(() => {
      const zs = [...document.querySelectorAll('#tafel-stempel tr')].slice(1);
      return zs.map(z => [z.cells[0].innerText, parseFloat(z.cells[1].innerText)]);
    });
    const sg = t.find(z => /Scharfetter/.test(z[0]));
    U(`Zweigstempel ${name}: Scharfetter-Gummel knotenweise exakt`,
      sg && sg[1] < 1e-11, sg ? sg[1].toExponential(2) : 'keine Tafel');
    if (sp === 43) {
      const ze = t.find(z => /zentral/.test(z[0]));
      U('Zweigstempel: die zentrale Differenz liegt dort weit daneben',
        ze && ze[1] > 0.05, ze ? ze[1].toExponential(2) : '—');
    }
  }

  // ---- Vergleich 2 und 3: Sprachproben ---------------------------------------------
  for (const [ix, nm] of [[2, 'Drucknetz'], [3, 'Impulsoperator']]) {
    await seite.click(`#wahl-vergleich button:nth-child(${ix})`);
    await seite.waitForTimeout(200);
    const wahlDa = await seite.evaluate(() => !document.getElementById('wahl-koerper').hidden);
    if (!wahlDa) {
      const grund = await seite.evaluate(() => document.getElementById('h-kein-koerper').innerText.trim());
      U(`${nm}: der Störkörper-Knopf fehlt UND die Seite sagt warum`,
        grund.length > 40, grund.split('\n')[0]);
      const z = await seite.evaluate(() => {
        const zs = [...document.querySelectorAll('#tafel-netz tr')].slice(1);
        return zs.map(r => r.cells[3].innerText.trim());
      });
      U(`${nm}: alle Proben RICHTIG`, z.length > 0 && z.every(u => u === 'RICHTIG'),
        `${z.filter(u => u === 'RICHTIG').length}/${z.length}`);
      continue;
    }
    const gesehen = [];
    for (const koerper of [1, 0]) {   // nth-child(1) ist „mit Störkörper"
      await seite.click(`#wahl-koerper button:nth-child(${koerper ? 1 : 2})`);
      await seite.waitForTimeout(700);
      const z = await seite.evaluate(() => {
        const zs = [...document.querySelectorAll('#tafel-netz tr')].slice(1);
        return { urteile: zs.map(r => r.cells[3].innerText.trim()),
                 werte: zs.map(r => r.cells[1].innerText.trim()),
                 zeilen: zs.map(r => r.cells[0].innerText.trim()).join(' | ') };
      });
      const alle = z.urteile.length > 0 && z.urteile.every(u => u === 'RICHTIG');
      U(`${nm}, ${koerper ? 'mit' : 'ohne'} Störkörper: alle Proben RICHTIG`,
        alle, `${z.urteile.filter(u => u === 'RICHTIG').length}/${z.urteile.length} · ${z.werte[0]}`);
      gesehen.push(z.zeilen + ' :: ' + z.werte.join(','));
    }
    // Die beiden Einstellungen MUESSEN verschiedene Proben ergeben -- taeten sie es
    // nicht, waere der Knopf wirkungslos. Genau so ist Befund B11 aufgefallen.
    U(`${nm}: „mit Störkörper" und „nacktes Netz" sind wirklich verschieden`,
      gesehen.length === 2 && gesehen[0] !== gesehen[1],
      gesehen.length === 2 && gesehen[0] !== gesehen[1] ? 'verschiedene Proben' : 'IDENTISCH');
  }

  // ---- Vergleich 4: der ganze Loeser, beide Wege ------------------------------------
  await seite.click('#wahl-vergleich button:nth-child(4)');
  await seite.waitForTimeout(200);
  await seite.evaluate(() => {
    document.getElementById('s-netz').value = '110,21';
    document.getElementById('s-tende').value = '4';      // 0,4 s
    zeigeTende();
  });
  await seite.click('#b-start');
  await seite.waitForFunction(() => !document.getElementById('b-start').disabled, { timeout: 180000 });
  const L = await seite.evaluate(() => ({
    du: lauf.du, dcD: lauf.dcD, t: lauf.ZA.t, k: lauf.ZA.k,
    cDA: lauf.ZA.cD, cDB: lauf.ZB.cD
  }));
  U('Ganzer Löser lief beide Wege durch', L.k > 100, `${L.k} Zeitschritte bis t = ${L.t.toFixed(3)} s`);
  U('beide Wege ergeben dasselbe Geschwindigkeitsfeld', L.du < 1e-10, L.du.toExponential(3) + ' relativ');
  U('beide Wege ergeben denselben Widerstandsbeiwert', L.dcD < 1e-8, L.dcD.toExponential(3));

  // ---- Zeitverlauf -----------------------------------------------------------------
  await seite.click('#r-zeitverlauf'); await seite.waitForTimeout(500);
  const hatZeit = await seite.evaluate(() => !!(window.LABORDATEN && window.LABORDATEN.zeit));
  if (hatZeit) {
    const toene = () => seite.evaluate(() => {
      const g = n => { const c = document.getElementById(n);
        const x = c.getContext('2d').getImageData(0, 0, c.width, c.height).data;
        const s = new Set();
        for (let k = 0; k < x.length; k += 1200) s.add(x[k] + ',' + x[k+1] + ',' + x[k+2]);
        return s.size; };
      return { z1: g('c-z1'), z2: g('c-z2'), z3: g('c-z3') };
    });
    for (const art of ['unabhaengig', 'wege', 'surrogat']) {
      await seite.evaluate(a => { waehleZeit(a);
        document.getElementById('s-lupe').value = '0';
        document.getElementById('s-zeit').value = String(Math.floor(+document.getElementById('s-zeit').max / 3));
        zeigeZeit(+document.getElementById('s-zeit').value); }, art);
      await seite.waitForTimeout(350);
      const z = await toene();
      U(`Zeitverlauf „${art}": die beiden oberen Schirme zeigen ein Feld`,
        z.z1 > 8 && z.z2 > 8, `Farbtöne ${z.z1}/${z.z2}`);
      // Bei GLEICHER Skala muss der Unterschied zweier WEGE verschwinden -- das ist
      // die Aussage. Beim SURROGAT dagegen soll er sichtbar sein: Sein Fehler ist
      // echt und wird nicht versteckt. Zwei Faelle, zwei Kriterien.
      // (Der Schirm ist nie ganz einfarbig: Der Zylinder liegt darauf.)
      if (art === 'unabhaengig') {
        // Farbtoene zu zaehlen war hier der falsche Massstab. Die Felder zweier
        // verschieden hergeleiteter Verfahren unterscheiden sich punktweise um rund
        // ein Viertel -- das ist gemessen und steht auf der Seite. Geprueft wird
        // deshalb die MESSZAHL, nicht der Augenschein: Die Differenz muss kleiner
        // sein als das Feld, und die Seite muss sagen, wie viel kleiner.
        const m = await seite.evaluate(() => window.LABORDATEN.zeit.galerkin.masse);
        U('Zeitverlauf „unabhaengig": die Differenz ist kleiner als das Feld',
          m.gesamt < 50, `Effektivwert ${m.gesamt.toFixed(0)} % des Feldes`);
        U('Zeitverlauf „unabhaengig": nahe am Körper größer als im Nachlauf',
          m.nahe > m.nachlauf,
          `Nachlauf ${m.nachlauf.toFixed(0)} %, nahe am Zylinder ${m.nahe.toFixed(0)} %`);
      } else if (art === 'wege')
        U('Zeitverlauf „wege": der Unterschied verschwindet bei gleicher Skala',
          z.z3 * 50 < z.z1, `${z.z3} gegen ${z.z1} Farbtöne`);
      else
        // Nicht „sichtbar" ist das Kriterium -- das war eine Erwartung aus der Zeit, als
        // das Surrogat 15,9 % Fehler hatte. Worauf es ankommt: Der Fehler muss KLEIN
        // sein und BEZIFFERT werden. Dass er bei gleicher Skala verschwindet, ist das
        // bessere Ergebnis, nicht das schlechtere.
        U('Zeitverlauf „surrogat": der Fehler ist kleiner als das Feld',
          z.z3 * 3 < z.z1, `${z.z3} gegen ${z.z1} Farbtöne`);
      // Und bei hoher Verstaerkung muss etwas erscheinen, sonst waere der Regler tot.
      // Der Verstaerkungsregler ist nur dort sinnvoll, wo die Differenz WEIT unter der
      // Skala liegt. Beim unabhaengigen Vergleich fuellt sie die Skala bereits aus --
      // dort waere „mehr Verstaerkung" ohne Bedeutung, und die Seite sagt das auch.
      if (art !== 'unabhaengig') {
        await seite.evaluate(() => { document.getElementById('s-lupe').value = '170'; zeigeZeit(null); });
        await seite.waitForTimeout(250);
        const z2 = await toene();
        U(`Zeitverlauf „${art}": der Verstärkungsregler wirkt`,
          z2.z3 > z.z3 + 5, `${z.z3} -> ${z2.z3} Farbtöne`);
      }
      await seite.evaluate(() => { document.getElementById('s-lupe').value = '0'; zeigeZeit(null); });
    }
    const gal = await seite.evaluate(() => {
      waehleZeit('unabhaengig');
      return { da: !!(window.LABORDATEN.zeit && window.LABORDATEN.zeit.galerkin),
               text: document.getElementById('h-zeit').innerText };
    });
    U('die Seite nennt die Beiwerte des unabhaengigen Wegs',
      gal.da && /3,2\d\d|3\.2\d\d/.test(gal.text), (gal.text.split('\n').pop() || '').slice(0, 70));
    const nennt = await seite.evaluate(() => {
      waehleZeit('surrogat');
      return document.getElementById('h-zeitfuss').innerText + ' ' +
             document.getElementById('h-zeit').innerText;
    });
    U('die Seite nennt den gemessenen Surrogatfehler in Zahlen',
      /\d+[,.]\d+\s*%/.test(nennt), (nennt.match(/[^\n]*%[^\n]*/) || ['keine Zahl'])[0].trim());
    const be = await seite.evaluate(() => window.LABORDATEN.zeit.bericht);
    U('der gemessene Surrogatfehler liegt unter 1 % der Schwankungsbreite',
      be.fehler_harmonisch < 0.01, `${(100*be.fehler_harmonisch).toFixed(3)} % ` +
      `(linearer Zeitoperator zum Vergleich: ${(100*be.fehler_dmd).toFixed(2)} %)`);
    U('die beiden Surrogate erben die Übereinstimmung der beiden Wege',
      be.unterschied_surrogat < 1e-9,
      `${be.unterschied_surrogat.toExponential(2)} gegen ${be.unterschied_echt.toExponential(2)} echt`);
    await seite.evaluate(() => waehleZeit('wege'));
    await seite.evaluate(() => { document.getElementById('s-zeit').value = '0'; zeigeZeit(0); });
    const f0 = await seite.evaluate(() => document.getElementById('w-zeit').innerText);
    await seite.evaluate(() => { const s = document.getElementById('s-zeit');
      s.value = String(Math.floor(+s.max / 2)); zeigeZeit(+s.value); });
    const f1 = await seite.evaluate(() => document.getElementById('w-zeit').innerText);
    U('der Zeitschieber bewegt den Verlauf', f0 !== f1, `${f0} -> ${f1}`);
  } else {
    U('Zeitdaten eingebaut', false, 'window.LABORDATEN.zeit fehlt');
  }

  // ---- Bilder zum Ansehen -----------------------------------------------------------
  for (const r of ['labor', 'zeitverlauf', 'wirbel', 'verfahren', 'befunde', 'anleitung']) {
    await seite.click('#r-' + r); await seite.waitForTimeout(500);
    await seite.screenshot({ path: path.join(BILD, `labor_${r}.png`), fullPage: false });
  }
  U('keine Konsolenfehler ueber den ganzen Durchlauf', fehler.length === 0,
    fehler.slice(0, 3).join(' | ') || 'keine');

  // ---- Galerkin: beide Wege rechnen WIRKLICH in der Seite (K15) --------------------
  {
    await seite.click('#r-labor'); await seite.waitForTimeout(200);
    await seite.click('#wahl-vergleich button:nth-child(5)');
    await seite.waitForTimeout(300);
    U('Galerkin-Vergleich: Einstellkasten geht auf',
      await seite.evaluate(() => !document.getElementById('k-galerkin').hidden), 'sichtbar');
    U('Triangle ist eingebettet, nicht nachgeladen',
      await seite.evaluate(() => (window.TRIANGLE_WASM || '').length > 100000 &&
                                 typeof window.TRIANGLE === 'object'),
      await seite.evaluate(() => `${Math.round((window.TRIANGLE_WASM||'').length/1024)} kB als base64`));
    // kurzer Lauf -- es geht um das RECHNEN, nicht um die Dauer
    await seite.evaluate(() => { document.getElementById('s-gtende').value = '5'; zeigeGTende(); });
    await seite.click('#b-gstart');
    let fertig = false;
    for (let i = 0; i < 200; i++) {
      await seite.waitForTimeout(1000);
      fertig = await seite.evaluate(() => !document.getElementById('b-gstart').disabled);
      if (fertig) break;
    }
    U('beide Wege sind durchgelaufen', fertig,
      await seite.evaluate(() => document.getElementById('b-gbahn').textContent));
    const z = await seite.evaluate(() => {
      const txt = document.getElementById('tafel-galerkin').innerText;
      const zahlen = [...txt.matchAll(/-?\d+\.\d+/g)].map(Number);
      const toene = id => { const c = document.getElementById(id);
        const d = c.getContext('2d').getImageData(0,0,c.width,c.height).data;
        const s = new Set(); for (let i=0;i<d.length;i+=4) s.add(d[i]+','+d[i+1]+','+d[i+2]);
        return s.size; };
      return { zahlen: zahlen, hinweis: document.getElementById('h-galerkin').textContent,
               fazit: document.getElementById('h-gfazit').textContent,
               netzbild: toene('c-gnetz'), gfeld: toene('c-gfeld'),
               nfeld: toene('c-nfeld'), bahn: toene('c-gbahn') };
    });
    U('das Netz wird gezeichnet (Dreiecke, Zylinderrand)', z.netzbild > 2,
      `${z.netzbild} Farbtöne`);
    /* Shewchuk bittet ausdruecklich: „If you use a mesh generated by Triangle in a
     * publication, please include an acknowledgment as well. And please spell Triangle
     * with a capital ``T''!" Das ist eine Bedingung des Urhebers, keine Hoeflichkeit --
     * und eine Bedingung, die unbemerkt aus einer Seite verschwinden kann. Deshalb ist
     * sie ein Pruefposten: SICHTBAR, mit dem Namen des Autors, und mit grossem T. */
    const dank = await seite.evaluate(() => {
      const d = document.querySelector('#f-galerkin .dank');
      if (!d) return { da: false };
      const r = d.getBoundingClientRect();
      return { da: true, sichtbar: r.height > 10 && r.width > 100,
               text: (d.innerText || '').replace(/\s+/g, ' ') };
    });
    U('die Danksagung an den Urheber des Vernetzers steht SICHTBAR in der Seite',
      dank.da && dank.sichtbar && /Triangle/.test(dank.text) &&
      /Shewchuk/.test(dank.text) && /Dank/.test(dank.text) &&
      /nicht verkauft|kommerzielle/.test(dank.text),
      dank.da ? dank.text.slice(0, 100) : 'kein Dankkasten gefunden');
    /* Befund 01.11.2026: Hier stand /\btriangle\b/ -- und das traf `triangle-wasm`.
     * Das ist aber ein PAKETNAME und gehoert klein geschrieben; die Bitte des Autors
     * gilt dem PROGRAMM. Die Probe meldete einen Fehler, den es nicht gab. Dateinamen
     * und Paketnamen werden deshalb vorher herausgenommen. */
    const ohneNamen = dank.da
      ? dank.text.replace(/triangle-wasm/g, '').replace(/triangle\.[a-z]+/g, '')
                 .replace(/labor\/triangle/g, '')
      : '';
    U('Triangle wird mit grossem T geschrieben (ausdrückliche Bitte des Autors)',
      dank.da && !/\btriangle\b/.test(ohneNamen),
      dank.da ? 'großes T überall, wo das Programm gemeint ist' : '—');
    U('beide Wirbelstärkefelder zeigen ein Feld, nicht eine Fläche',
      z.gfeld > 200 && z.nfeld > 200, `Galerkin ${z.gfeld}, Netzwerk ${z.nfeld} Farbtöne`);
    /* Befund 01.11.2026: Hier stand nur `z.bahn > 50`, also eine Zahl von Farbtoenen.
     * Die war erfuellt, obwohl im Diagramm nur EINE der beiden Kurven zu sehen war --
     * die Legende versprach zwei. Eine Farbzahl sagt nicht, WELCHE Farben da sind.
     * Jetzt wird nach den beiden Kurvenfarben selbst gesucht: Blau fuer Galerkin,
     * Rot fuer das Netzwerk. Gefunden hat den Fehler das ANSEHEN des Bildes, nicht
     * die Probe (Hausregel „Bild ansehen, nicht nur Pruefzahl"). */
    const kurven = await seite.evaluate(() => {
      const c = document.getElementById('c-gbahn');
      const d = c.getContext('2d').getImageData(0, 0, c.width, c.height).data;
      let blau = 0, rot = 0;
      for (let i = 0; i < d.length; i += 4) {
        const r = d[i], g = d[i+1], b = d[i+2];
        if (b > r + 30 && b > 100 && r < 140) blau++;
        if (r > b + 40 && r > 120 && b < 120) rot++;
      }
      return { blau: blau, rot: rot };
    });
    U('BEIDE Kurven sind im Zeitdiagramm gezeichnet, nicht nur eine',
      kurven.blau > 200 && kurven.rot > 200,
      `blau (Galerkin) ${kurven.blau} Punkte, rot (Netzwerk) ${kurven.rot} Punkte`);
    const cD_G = z.zahlen[0], cD_N = z.zahlen[1];
    U('Galerkin liefert einen plausiblen Widerstandsbeiwert',
      isFinite(cD_G) && cD_G > 1.5 && cD_G < 6, `c_D = ${cD_G}`);
    U('Netzwerk liefert einen plausiblen Widerstandsbeiwert',
      isFinite(cD_N) && cD_N > 1.5 && cD_N < 6, `c_D = ${cD_N}`);
    /* Die beiden duerfen sich unterscheiden -- sie sind verschieden hergeleitet und
     * rechnen auf grobem Netz. Was NICHT sein darf: dass sie identisch sind. Dann
     * waere einer von beiden nicht gerechnet worden, sondern abgeschrieben. */
    U('die beiden Wege sind wirklich zwei Rechnungen, keine Kopie',
      Math.abs(cD_G - cD_N) > 1e-6, `Unterschied ${Math.abs(cD_G-cD_N).toFixed(4)}`);
    U('die Seite nennt den Knotenabstand beider Wege (faire Gegenüberstellung)',
      /Knoten alle/.test(z.hinweis) && /gleicher Knotenabstand/.test(z.hinweis),
      z.hinweis.split('\n')[0].slice(0, 80));
    U('die Seite sagt, dass die Browserwerte nicht die Benchmarkwerte sind',
      /3,2318/.test(z.fazit) && /grobem Netz|groben Netz/.test(z.fazit),
      z.fazit.split('\n')[2] || z.fazit.slice(0, 80));
  }

  await browser.close();
  console.log('='.repeat(96));
  console.log(' KNOPFPROBE CFD-LABOR -- ' + path.basename(DATEI));
  console.log('='.repeat(96));
  let n = 0;
  for (const [t, ok, b] of urteile) { if (ok) n++; console.log(`  [${ok ? 'RICHTIG' : 'FALSCH '}] ${t}\n              ${b}`); }
  console.log(`\n  == ${n}/${urteile.length}`);
  process.exit(n === urteile.length ? 0 : 1);
})();
