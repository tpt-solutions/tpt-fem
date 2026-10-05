/* Steady heat conduction on a unit cube: u = 0 at x = 0, u = 1 at x = 1, no
 * source, so the exact solution is u = x. Prints the max error.
 *
 *   cargo build --release
 *   gcc examples/poisson.c -Iinclude -Ltarget/release -ltpt_fem_capi -o poisson
 *   (Linux: LD_LIBRARY_PATH=target/release ./poisson)
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "tpt_fem.h"

static size_t *nodes_on_x(const TptMesh *mesh, double x, size_t *count) {
    size_t *ids;
    tpt_mesh_nodes_on_plane(mesh, 0, x, 1e-9, NULL, 0, count); /* size query */
    ids = malloc(*count * sizeof *ids);
    tpt_mesh_nodes_on_plane(mesh, 0, x, 1e-9, ids, *count, count);
    return ids;
}

int main(void) {
    const double lo[3] = {0, 0, 0}, hi[3] = {1, 1, 1};
    const size_t n[3] = {4, 4, 4};
    TptMesh *mesh = NULL;
    TptField *u = NULL;
    size_t n0, n1, nb, i, nn;
    size_t *left, *right;
    TptNodeBc *bcs;
    double max_err = 0.0;

    printf("tpt-fem %s\n", tpt_version());
    if (tpt_mesh_box(lo, hi, n, &mesh) != TPT_OK) goto fail;

    left = nodes_on_x(mesh, 0.0, &n0);
    right = nodes_on_x(mesh, 1.0, &n1);
    nb = n0 + n1;
    bcs = malloc(nb * sizeof *bcs);
    for (i = 0; i < n0; ++i) { bcs[i].node = left[i];       bcs[i].value = 0.0; }
    for (i = 0; i < n1; ++i) { bcs[n0 + i].node = right[i]; bcs[n0 + i].value = 1.0; }

    if (tpt_solve_poisson(mesh, 1.0, 2, 0.0, NULL, NULL, bcs, nb, &u) != TPT_OK) goto fail;

    nn = tpt_mesh_node_count(mesh);
    for (i = 0; i < nn; ++i) {
        double c[3], e;
        tpt_mesh_coords(mesh, i, c);
        e = fabs(tpt_field_data(u)[i] - c[0]);
        if (e > max_err) max_err = e;
    }
    printf("%zu nodes, max |u - x| = %.3e\n", nn, max_err);

    free(left); free(right); free(bcs);
    tpt_field_free(u);
    tpt_mesh_free(mesh);
    return max_err < 1e-8 ? 0 : 1;

fail:
    fprintf(stderr, "error: %s\n", tpt_last_error_message());
    return 2;
}
