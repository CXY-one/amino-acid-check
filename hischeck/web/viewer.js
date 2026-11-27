"use strict";

window.hisViewer = (() => {
    const panel = document.getElementById("scene");
    const empty = document.getElementById("empty");
    const focus = document.getElementById("focus");
    const neighbors = document.getElementById("neighbors");
    const evidence = document.getElementById("evidence");
    const nearby = document.getElementById("nearby");
    const measure = document.getElementById("measure");
    const clearMeasurements = document.getElementById("clear-measurements");
    const legend = document.getElementById("legend");
    const caption = document.getElementById("caption");
    const measureColor = "#245ac4";
    const coordinationColor = "#8854a3";
    let text = {};
    const tr = (key, values = {}) => (text[`viewer.${key}`] || key).replace(
        /\{(\w+)\}/g, (placeholder, name) => values[name] ?? placeholder);
    const elements = { C: "#a9d46c", N: "#258cff", O: "#cb4450", H: "#fffefe", D: "#fffefe",
        S: "#c59729", Se: "#c59729", P: "#b86f27" };
    const neighborElements = { ...elements, C: "#bed9a0", N: "#4b9bec", O: "#bd6671", H: "#edf0f3", D: "#edf0f3" };
    const neighborStyle = { stick: { radius: 0.125, colorscheme: { map: neighborElements, prop: "elem" } } };
    const activeStyle = { stick: { radius: 0.13, colorscheme: { map: elements, prop: "elem" } } };
    let viewer;
    let model;
    let structure;
    let atomBySerial;
    let atomByKey;
    let homeView;
    let evidenceAtoms;
    let measuring = false;
    let firstAtom = null;
    let measurements = [];
    let captionText = () => tr("hint");

    const chainOf = (atom) => atom.chain?.trim() || "_";
    const residueOf = (atom) => `${chainOf(atom)}:${atom.resi}${atom.icode?.trim() || ""}`;
    const atomKey = (atom) => JSON.stringify([chainOf(atom), atom.resi, atom.icode?.trim() || "", atom.atom]);
    const atomName = (atom) => `${residueOf(atom)} ${atom.resn} ${["H", "D"].includes(atom.elem) ? atom.elem : atom.atom}`;
    const point = (atom) => ({ x: atom.x, y: atom.y, z: atom.z });
    const distance = (first, second) => Math.hypot(first.x - second.x, first.y - second.y, first.z - second.z);
    const pairKey = (first, second) => first < second ? `${first}:${second}` : `${second}:${first}`;

    function showCaption(format) {
        captionText = format;
        caption.textContent = captionText();
    }

    function setMeasuring(enabled) {
        measuring = enabled;
        firstAtom = null;
        measure.setAttribute("aria-pressed", String(enabled));
        panel.classList.toggle("measuring", enabled);
        showCaption(() => tr(enabled ? "first_atom" : "hint"));
    }

    function pickAtom(atom) {
        if (!measuring) {
            showCaption(() => `${atomName(atom)} · ${atom.elem}`);
            return;
        }
        const key = atomKey(atom);
        if (firstAtom === null || firstAtom === key) {
            firstAtom = key;
            showCaption(() => tr("next_atom", { atom: atomName(atom) }));
        } else {
            const from = atomByKey.get(firstAtom);
            if (!measurements.some(([a, b]) => (a === firstAtom && b === key) || (a === key && b === firstAtom))) {
                measurements.push([firstAtom, key]);
            }
            firstAtom = null;
            showCaption(() => `${atomName(from)} ↔ ${atomName(atom)} · ${distance(from, atom).toFixed(2)} Å`
                + ` · ${tr("measured", { count: measurements.length })}`);
        }
        render();
    }

    function setLocale(locale) {
        text = locale.text;
        document.documentElement.lang = locale.language;
        for (const element of document.querySelectorAll("[data-text]")) {
            element.textContent = tr(element.dataset.text);
        }
        for (const element of document.querySelectorAll("[data-tip]")) {
            element.title = tr(element.dataset.tip);
        }
        panel.setAttribute("aria-label", tr(panel.dataset.aria));
        caption.textContent = captionText();
        render();
    }

    function initialize(locale) {
        if (locale) setLocale(locale);
        if (viewer) return "";
        try {
            viewer = $3Dmol.createViewer(panel, { backgroundColor: "#f2f6f9", antialias: true });
            viewer.setProjection("orthographic");
            viewer.render();
            caption.hidden = true;
            return "";
        } catch (error) {
            viewer = undefined;
            return tr("init_error", { error: error.message || error });
        }
    }

    function clear() {
        if (viewer) {
            viewer.removeAllLabels();
            viewer.removeAllShapes();
            viewer.removeAllModels();
            viewer.render();
        }
        structure = model = atomBySerial = atomByKey = homeView = evidenceAtoms = undefined;
        measurements = [];
        setMeasuring(false);
        focus.disabled = neighbors.disabled = evidence.disabled = nearby.disabled = measure.disabled = true;
        clearMeasurements.disabled = true;
        caption.hidden = true;
        empty.hidden = false;
        legend.replaceChildren();
    }

    function orientation([x, y]) {
        // The rows give the horizontal, vertical and normal axes of the target group.
        const [xx, xy, xz] = x;
        const [yx, yy, yz] = y;
        const normal = new $3Dmol.Vector3().crossVectors(new $3Dmol.Vector3(...x), new $3Dmol.Vector3(...y));
        const { x: zx, y: zy, z: zz } = normal;
        const trace = xx + yy + zz;
        // Choose the largest quaternion component to avoid cancellation near 180 degrees.
        if (trace > 0) {
            const scale = 2 * Math.sqrt(1 + trace);
            return [(zy - yz) / scale, (xz - zx) / scale, (yx - xy) / scale, scale / 4];
        }
        if (xx >= yy && xx >= zz) {
            const scale = 2 * Math.sqrt(1 + xx - yy - zz);
            return [scale / 4, (xy + yx) / scale, (xz + zx) / scale, (zy - yz) / scale];
        }
        if (yy >= zz) {
            const scale = 2 * Math.sqrt(1 - xx + yy - zz);
            return [(xy + yx) / scale, scale / 4, (yz + zy) / scale, (xz - zx) / scale];
        }
        const scale = 2 * Math.sqrt(1 - xx - yy + zz);
        return [(xz + zx) / scale, (yz + zy) / scale, scale / 4, (yx - xy) / scale];
    }

    function viewAroundTarget() {
        const selection = [...structure.focusAtoms, ...structure.hydrogens];
        const visible = [...selection, ...evidenceAtoms];
        viewer.zoomTo({ serial: visible.length ? visible : structure.target });
        const view = viewer.getView();
        const [cx, cy, cz] = structure.center;
        viewer.setView([-cx, -cy, -cz, view[3], ...orientation(structure.axes)]);
        viewer.rotate(13, "y");
        viewer.zoom(0.80);
        return viewer.getView();
    }

    function addLabel(text, position, color, options = {}) {
        viewer.addLabel(text, {
            position, fontSize: 13, fontColor: color, backgroundColor: "white", backgroundOpacity: 0.90,
            borderThickness: 0, inFront: true, ...options,
        });
    }

    function drawDistance(from, to, value, color, manual = false) {
        viewer.addCylinder({ start: point(from), end: point(to), color, radius: manual ? 0.05 : 0.039,
            dashed: true, dashLength: 0.15, gapLength: 0.09 });
        const fraction = manual ? 0.50 : 0.64;
        const at = (axis) => from[axis] + fraction * (to[axis] - from[axis]);
        addLabel(`${value.toFixed(2)} Å`, { x: at("x"), y: at("y"), z: at("z") }, color,
            { fontSize: 12, alignment: "center", screenOffset: { x: 0, y: manual ? -32 : 24 } });
    }

    function addLegend(kind, color) {
        const item = document.createElement("span");
        const sample = document.createElement("i");
        sample.className = "line";
        sample.style.borderColor = color;
        item.append(sample, tr(`contact.${kind}`));
        legend.append(item);
    }

    function render() {
        if (!model) return;
        viewer.removeAllLabels();
        viewer.removeAllShapes();
        model.setStyle({}, neighbors.checked ? neighborStyle : {});
        if (neighbors.checked) {
            model.setStyle({ serial: structure.waters }, { sphere: { radius: 0.18, color: "#d79aa2" } });
        }
        const metals = new Set(Object.keys(structure.metals).map(Number));
        const contactNonmetals = evidenceAtoms.filter(serial => !metals.has(serial));
        model.setStyle({ serial: [...structure.target, ...contactNonmetals] }, activeStyle);
        model.setStyle({ serial: [...structure.focusAtoms, ...structure.hydrogens] }, {
            stick: { radius: 0.16, colorscheme: { prop: "elem", map: elements } },
            sphere: { scale: 0.20, colorscheme: { prop: "elem", map: elements } },
        });
        for (const serial of evidenceAtoms) {
            const atom = atomBySerial.get(serial);
            if (atom && !metals.has(serial) && !structure.target.includes(serial)) {
                model.setStyle({ serial }, { sphere: { radius: 0.26, color: elements[atom.elem] || "#50778b" } }, true);
            }
        }
        for (const serial of metals) {
            const atom = atomBySerial.get(serial);
            if (!atom || !(neighbors.checked || evidenceAtoms.includes(serial) ||
                (nearby.checked && structure.nearbyAtoms.includes(serial)))) continue;
            model.setStyle({ serial }, { sphere: { radius: 0.42, color: "#8660a2" } });
            addLabel(structure.metals[serial], atom, "#693b80", { alignment: "center", screenOffset: { x: 0, y: -18 } });
        }
        if (nearby.checked) {
            for (const serial of structure.nearbyAtoms) {
                const atom = atomBySerial.get(serial);
                if (!atom || metals.has(serial)) continue;
                if (!evidenceAtoms.includes(serial)) {
                    model.setStyle({ serial }, { sphere: { radius: 0.23,
                        color: elements[atom.elem] || "#50778b" } }, true);
                }
                addLabel(`${residueOf(atom)} ${atom.atom}`, atom, "#3d5366",
                    { fontSize: 10, alignment: "center", screenOffset: { x: 0, y: 16 } });
            }
        }
        for (const serial of structure.labels) {
            const atom = atomBySerial.get(serial);
            if (atom) addLabel(atom.atom, atom, atom.elem === "N" ? "#1b59a0" : "#254758");
        }
        legend.replaceChildren();
        const pairs = measurements.map(([first, second]) => [atomByKey.get(first), atomByKey.get(second)]);
        const manualPairs = new Set(pairs.map(([from, to]) => pairKey(from.serial, to.serial)));
        if (evidence.checked) {
            const shown = new Set();
            for (const link of structure.links) {
                const from = atomBySerial.get(link.from);
                const to = atomBySerial.get(link.to);
                if (!from || !to || manualPairs.has(pairKey(from.serial, to.serial))) continue;
                drawDistance(from, to, link.distance, link.color);
                if (link.kind === "hydrogen" || shown.has(link.kind)) continue;
                shown.add(link.kind);
                addLegend(link.kind, link.color);
            }
            if (neighbors.checked) {
                let hasCoordination = false;
                for (const link of structure.coordination) {
                    const from = atomBySerial.get(link.from);
                    const to = atomBySerial.get(link.to);
                    if (!from || !to || manualPairs.has(pairKey(from.serial, to.serial))) continue;
                    viewer.addCylinder({ start: point(from), end: point(to), color: coordinationColor,
                        radius: 0.047, dashed: true, dashLength: 0.15, gapLength: 0.09 });
                    hasCoordination = true;
                }
                if (hasCoordination && !shown.has("metal")) addLegend("metal", coordinationColor);
            }
        }
        for (const [from, to] of pairs) drawDistance(from, to, distance(from, to), measureColor, true);
        if (pairs.length) addLegend("manual", measureColor);
        const endpoints = new Set(pairs.flat());
        if (firstAtom !== null) endpoints.add(atomByKey.get(firstAtom));
        for (const atom of endpoints) {
            viewer.addSphere({ center: point(atom), radius: 0.31, color: measureColor, opacity: 0.42 });
        }
        clearMeasurements.disabled = !pairs.length && firstAtom === null;
        viewer.render();
    }

    function setBond(atom, other, order) {
        const index = atom.bonds.indexOf(other.index);
        if (index < 0) {
            atom.bonds.push(other.index);
            atom.bondOrder.push(order);
        } else {
            atom.bondOrder[index] = order;
        }
    }

    function load(scene) {
        const error = initialize();
        if (error) return error;
        try {
            const sameSite = structure?.site === scene.site;
            const previous = sameSite ? viewer.getView() : null;
            const defaultView = sameSite && homeView;
            const saved = sameSite ? measurements : [];
            const wasMeasuring = sameSite && measuring;
            clear();
            model = viewer.addModel(scene.pdb, "pdb", { keepH: true, noSecondaryStructure: true });
            const atoms = model.selectedAtoms({});
            atomBySerial = new Map(atoms.map(atom => [atom.serial, atom]));
            atomByKey = new Map(atoms.map(atom => [atomKey(atom), atom]));
            if (scene.target.some(serial => !atomBySerial.has(serial))) throw new Error(tr("missing_target"));
            const metalIndices = new Set(atoms.filter(atom => Object.hasOwn(scene.metals, atom.serial)).map(atom => atom.index));
            for (const atom of atoms) {
                for (let i = atom.bonds.length - 1; i >= 0; i--) {
                    if (metalIndices.has(atom.index) || metalIndices.has(atom.bonds[i])) {
                        atom.bonds.splice(i, 1);
                        atom.bondOrder.splice(i, 1);
                    }
                }
            }
            for (const [from, to, order] of scene.bonds) {
                if (Object.hasOwn(scene.metals, from) || Object.hasOwn(scene.metals, to)) continue;
                const first = atomBySerial.get(from);
                const second = atomBySerial.get(to);
                if (first && second) {
                    setBond(first, second, order);
                    setBond(second, first, order);
                }
            }
            structure = scene;
            evidenceAtoms = [...new Set(scene.links.map(link => link.to))];
            measurements = saved.filter(([a, b]) => atomByKey.has(a) && atomByKey.has(b));
            model.setClickable({}, true, pickAtom);
            setMeasuring(wasMeasuring);
            render();
            homeView = defaultView || viewAroundTarget();
            if (previous) viewer.setView(previous);
            focus.disabled = neighbors.disabled = nearby.disabled = measure.disabled = false;
            evidence.disabled = scene.links.length === 0 && scene.coordination.length === 0;
            empty.hidden = true;
            caption.hidden = false;
            return "";
        } catch (error) {
            clear();
            return tr("scene_error", { error: error.message || error });
        }
    }

    focus.addEventListener("click", () => { if (homeView) viewer.setView(homeView); });
    neighbors.addEventListener("change", render);
    evidence.addEventListener("change", render);
    nearby.addEventListener("change", render);
    measure.addEventListener("click", () => {
        setMeasuring(!measuring);
        render();
    });
    clearMeasurements.addEventListener("click", () => {
        measurements = [];
        setMeasuring(measuring);
        render();
    });
    return { initialize, setLocale, clear, load };
})();
