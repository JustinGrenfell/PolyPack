// Imported by the tests that read polytopes with getPolytope: loads the entries kept on demand (the tables and 8:421).
import { LAZY_CODES, loadPolytope } from "@polypack/ndim";

await Promise.all(LAZY_CODES.map((c) => loadPolytope(c)));
