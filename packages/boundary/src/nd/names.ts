/**
 * The simplex, orthoplex and demicube as containers of the boundary space for d = 5..10 ("4D on 5D" to "9D on 10D",
 * docs/spaces/BOUNDARY.md section 4C): names, sizes and the helper that the registries share with the d-cubes.
 *
 * A container name is `<d><family>`: "5simplex" .. "10simplex", "5orthoplex" .. "10orthoplex", "5demicube" .. "10demicube". The d-cubes of
 * section 4B keep their names ("5cube" .. "10cube"); with these, every dimension d = 5..10 has four containers.
 */

export const ND_FAMILIES = ["simplex", "orthoplex", "demicube"] as const;
export type NDFamily = (typeof ND_FAMILIES)[number];
export const ND_DIMS = [5, 6, 7, 8, 9, 10] as const;
export type NDDim = (typeof ND_DIMS)[number];
export type NDName = `${NDDim}${NDFamily}`;

/** The 18 container names, simplex first, then orthoplex, then demicube, each by d. */
export const ND_NAMES: readonly NDName[] = ND_FAMILIES.flatMap((f) => ND_DIMS.map((d) => `${d}${f}` as NDName));

const NAME_RE = /^(5|6|7|8|9|10)(simplex|orthoplex|demicube)$/;

/** The family and ambient dimension of a name, or null when it is none of the 18. */
export function ndParse(name: unknown): { family: NDFamily; d: NDDim } | null {
  if (typeof name !== "string") return null;
  const m = NAME_RE.exec(name);
  return m ? { family: m[2] as NDFamily, d: Number(m[1]) as NDDim } : null;
}

export function isNDName(x: unknown): x is NDName {
  return ndParse(x) !== null;
}

/** The ambient dimension d of an nD container name, or 0. */
export function ndDim(name: unknown): number {
  return ndParse(name)?.d ?? 0;
}

/** Number of cells (facets) of the container. */
export function ndCellCount(family: NDFamily, d: number): number {
  if (family === "simplex") return d + 1;
  if (family === "orthoplex") return 2 ** d;
  return 2 * d + 2 ** (d - 1);
}

/** Pieces: balls only (polytope pieces on these boundaries are a later phase). */
export const ND_PIECES = ["ball"] as const;

/**
 * The most pieces of each container's problem set. P = 1..24 for every container: the work caps (ND_MAX_WORK) refuse a
 * check that needs more, they do not shorten the range.
 */
export const ND_MAX_P = 24;

/** The titles of the containers (the problem picker's). */
export const ND_TITLES: Record<NDName, string> = Object.fromEntries(
  ND_NAMES.map((n) => {
    const { family, d } = ndParse(n)!;
    return [n, `${d}-${family}`];
  }),
) as Record<NDName, string>;
