// run this with Node 22.19 or higher, in a directory with p-limit installed, store as a .ts-file.
// command:  node index.ts <input-directory> <output-directory>

import pLimit from "p-limit";
import fs from "node:fs";
import path from "node:path";
import { exec as execWithCb } from "node:child_process";
import { promisify } from "node:util";
const exec = promisify(execWithCb);

const input = process.argv[2];
const output = process.argv[3];

if (!input || !output) {
  console.log("usage: node index.ts <input-dir> <output-dir>");
  process.exit();
}

const files = fs.readdirSync(input);

const limit = pLimit(8);

const freesasaTasks = files.map((file) =>
  limit(() => {
    const inputFile = path.join(input, file);
    const outputFile = path.join(
      output,
      file.replace(".pdb", ".freesasa.json"),
    );
    return exec(
      `freesasa --shrake-rupley -n 100 -t 2 --output ${outputFile} --format json ${inputFile}`,
    );
  }),
);

const rustsasaTasks = files.map((file) =>
  limit(() => {
    const inputFile = path.join(input, file);
    const outputFile = path.join(
      output,
      file.replace(".pdb", ".rustsasa.json"),
    );
    return exec(
      `~/.cargo/bin/rust-sasa ${inputFile} ${outputFile} --format=json`,
    );
  }),
);

const startFreesasa = performance.now();
await Promise.all(freesasaTasks);
console.log(
  "Freesasa parallel processes ",
  ((performance.now() - startFreesasa) / 1000).toFixed(2),
  " seconds",
);

const startRust = performance.now();
await Promise.all(rustsasaTasks);
console.log(
  "Rustsasa parallel processes ",
  ((performance.now() - startRust) / 1000).toFixed(2),
  " seconds",
);

const startRust2 = performance.now();
await exec(`~/.cargo/bin/rust-sasa ${input} ${output} --format=json`);
console.log(
  "Rustsasa single batch process ",
  ((performance.now() - startRust2) / 1000).toFixed(2),
  " seconds",
);
