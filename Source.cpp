#include <igraph.h>
#include <string>
#include "../../FileIO.h"
#include <vector>

void GRAPH2iGRAPH(igraph_t &graph, GRAPH &grph);
void save(const std::string dir, igraph_vector_t &membership, igraph_matrix_t &memberships, igraph_vector_t modularity);

int main() {
	
	GRAPH grph;

	std::string search_dir = "D:/Projects/Thesis/iGraphGenerator/iGraphGenerator/Watts_Strogatz2/N_640/K_3/P_0.007000/";
	std::string curr_dir;

	std::vector< std::string > graph_paths;

	
	std::string tmp_str;

	
	searchDirs(search_dir, ".mat", graph_paths);


	printf("Number of Graphs = %d\n", graph_paths.size());


	
	for (size_t g = 0; g < graph_paths.size(); g++) {
		curr_dir = fileDirectory(graph_paths[g].c_str());

		tmp_str = curr_dir.substr(search_dir.size());
		/*
		if (g % 1000 == 0) {
			printf("\tGRAPH %d / %d: %s\n", g + 1, graph_paths.size(), tmp_str.c_str());
		}
		*/

		//printf("Loading Graph\n");
		load_graph(graph_paths[g], grph);
		const size_t N = grph.NoV;

		igraph_t graph;
		igraph_empty(&graph, N, false);

		//Community Identification
		igraph_vector_t membership;
		igraph_matrix_t memberships;
		igraph_vector_t modularity;

		igraph_vector_init(&membership, N);
		igraph_matrix_init(&memberships, N, N);
		igraph_vector_init(&modularity, N);

		igraph_vector_t nbrlist;
		igraph_vector_init(&nbrlist, 0);

		GRAPH2iGRAPH(graph, grph);

		//Community Identification
		igraph_community_multilevel(&graph, nullptr, &membership, &memberships, &modularity);

		save(curr_dir, membership, memberships, modularity);

		COMMUNITY_DATA com;
		load_community(curr_dir, com);

		igraph_vector_destroy(&nbrlist);
		igraph_destroy(&graph);

		//Community Identification
		igraph_vector_destroy(&membership);
		igraph_matrix_destroy(&memberships);
		igraph_vector_destroy(&modularity);
	}
	
	printf("Communities Identification Complete.\n");

	return 0;
}

void GRAPH2iGRAPH(igraph_t &graph, GRAPH &grph) {
	igraph_destroy(&graph);

	igraph_integer_t j;

	igraph_adjlist_t adjlist;
	igraph_adjlist_init_empty(&adjlist, grph.NoV);

	igraph_vector_int_t *nbrs;

	for (size_t i = 0; i < grph.NoV; i++) {
		nbrs = igraph_adjlist_get(&adjlist, i);
		igraph_vector_int_resize(nbrs, grph.deglist[i]);

		for (size_t k = 0; k < grph.deglist[i]; k++) {
			j = (igraph_integer_t)grph.adjlist[i*grph.max_deg + k];
			VECTOR(*nbrs)[k] = j;
		}
	}

	igraph_adjlist(&graph, &adjlist, IGRAPH_ALL, true);

	igraph_adjlist_destroy(&adjlist);
}

void save(const std::string dir, igraph_vector_t &membership, igraph_matrix_t &memberships, igraph_vector_t modularity) {
	const std::uint32_t N = (std::uint32_t) igraph_vector_size(&membership);

	const std::uint32_t W = (std::uint32_t) igraph_matrix_ncol(&memberships);
	const std::uint32_t H = (std::uint32_t) igraph_matrix_nrow(&memberships);

	COMMUNITIES coms(N, W, H);

	//printf("N = %d | W = %d | H = %d\n", N, W, H);
	

	for (size_t i = 0; i < N; i++) {
		coms.membership[i] = (std::uint16_t) igraph_vector_e(&membership, i);
	}

	//printf("Number of Communities = %d\n", *std::max_element(coms.membership, coms.membership + N));

	for (size_t i = 0; i < H; i++) {
		coms.modularities[i] = (double)igraph_vector_e(&modularity, i);
		for (size_t j = 0; j < W; j++) {
			coms.memberships[i*W + j] = (std::uint16_t) MATRIX(memberships, i, j);
		}
	}

	save_community(dir, coms);
}