#include <unistd.h>
#include <cstdio>
#include <vector>
#include <string>
#include <map>

#include "scoup.h"

using namespace std;

int main_scoup(int argc, char* argv[]){
	int K, max_ite1, max_ite2;
	double alpha_min, alpha_max, t_min, t_max, sigma_squared_min, thresh;
	K = 1;
	max_ite1 = 1000;
	max_ite2 = 10000;
	alpha_min = 0.1;
	alpha_max = 100;
	t_min = 0.001;
	t_max = 2.0;
	sigma_squared_min = 0.1;
	thresh = 0.01;

	int ch;
	extern char *optarg;
	extern int optind, opterr;
	while((ch = getopt(argc, argv, "k:m:M:a:A:t:T:s:e:")) != -1){
		switch(ch){
		case 'k':
			K = atoi(optarg);
			break;
		case 'm':
			max_ite1 = atoi(optarg);
			break;
		case 'M':
			max_ite2 = atoi(optarg);
			break;
		case 'a':
			alpha_min = atof(optarg);
			break;
		case 'A':
			alpha_max = atof(optarg);
			break;
		case 't':
			t_min = atof(optarg);
			break;
		case 'T':
			t_max = atof(optarg);
			break;
		case 's':
			sigma_squared_min = atof(optarg);
			break;
		case 'e':
			thresh = atof(optarg);
			break;
		case ':':
			printf("invalid option\n");
			return 1;
		case '?':
			printf("invalid option\n");
			return 1;
		}
	}

	FILE *fin_expression, *fin_init, *fin_time, *fout_gene_para, *fout_cell_para, *fout_ll;
	if((fin_expression=fopen(argv[optind], "r")) == NULL){
		printf("cannot open expression data\n");
		return 1;
	}
	if((fin_init=fopen(argv[optind+1], "r")) == NULL){
		printf("cannot open initial distribution data\n");
		return 1;
	}
	if((fin_time=fopen(argv[optind+2], "r")) == NULL){
		printf("cannot open initial pseudo-time data\n");
		return 1;
	}
	if((fout_gene_para=fopen(argv[optind+3], "w")) == NULL){
		printf("cannot open Output_file1 (parameters related to gene and lineage)\n");
	}
	if((fout_cell_para=fopen(argv[optind+4], "w")) == NULL){
		printf("cannot open Output_file2 (parameters related to cell)\n");
	}
	if((fout_ll=fopen(argv[optind+5], "w")) == NULL){
		printf("cannot open Output_file3 (log-likelihood)\n");
	}

	int gene_num = atoi(argv[optind+6]);
	int cell_num = atoi(argv[optind+7]);

	Continuous_OU_process OU(gene_num, cell_num, K, max_ite1, max_ite2, alpha_min, alpha_max, t_min, t_max, sigma_squared_min, thresh);

	OU.Set_expression(fin_expression);
	if(OU.Set_initial_parameter(fin_init) == 1){
		return 1;
	}
	if(OU.Set_time(fin_time) == 1){
		return 1;
	}
	OU.Set_initial_parameter2();

	//initialization with K=1 model
	for(int i=0; i<gene_num; i++){
		OU.EM_for_null_model(i);
	}
	//initialization with random responsibility
	OU.Init_EM();
	//optimization
	OU.EM();

	//output
	OU.Print_cell_parameter(fout_cell_para);
	OU.Print_gene_parameter(fout_gene_para);
	OU.Print_ll(fout_ll);

	return 0;
}

int main_scoup_resume(int argc, char* argv[]){
  int K, max_ite1, max_ite2;
  double alpha_min, alpha_max, t_min, t_max, sigma_squared_min, thresh;
  K = 1;
  max_ite1 = 100;
  max_ite2 = 100;
  alpha_min = 0.1;
  alpha_max = 100;
  t_min = 0.001;
  t_max = 2.0;
  sigma_squared_min = 0.1;
  thresh = 0.01;

  int ch;
  extern char *optarg;
  extern int optind, opterr;
  while((ch = getopt(argc, argv, "k:m:M:a:A:t:T:s:e:")) != -1){
    switch(ch){
    case 'k':
      K = atoi(optarg);
      break;
    case 'm':
      max_ite1 = atoi(optarg);
      break;
    case 'M':
      max_ite2 = atoi(optarg);
      break;
    case 'a':
      alpha_min = atof(optarg);
      break;
    case 'A':
      alpha_max = atof(optarg);
      break;
    case 't':
      t_min = atof(optarg);
      break;
    case 'T':
      t_max = atof(optarg);
      break;
    case 's':
      sigma_squared_min = atof(optarg);
      break;
    case 'e':
      thresh = atof(optarg);
      break;
    case ':':
      printf("invalid option\n");
      return 1;
    case '?':
      printf("invalid option\n");
      return 1;
    }
  }

  FILE *fin_expression, *fin_init, *fin_opt_gene_para, *fin_opt_cell_para, *fout_gene_para, *fout_cell_para, *fout_ll;
  if((fin_expression=fopen(argv[optind], "r")) == NULL){
    printf("cannot open expression data\n");
    return 1;
  }
  if((fin_init=fopen(argv[optind+1], "r")) == NULL){
    printf("cannot open initial distribution data\n");
    return 1;
  }
  if((fin_opt_gene_para=fopen(argv[optind+2], "r")) == NULL){
    printf("cannot open semi optimized gene and lineage parameters\n");
    return 1;
  }
  if((fin_opt_cell_para=fopen(argv[optind+3], "r")) == NULL){
    printf("cannot open semi optimized cell parameters\n");
    return 1;
  }
  if((fout_gene_para=fopen(argv[optind+4], "w")) == NULL){
    printf("cannot open Output_file1 (parameters related to gene and lineage)\n");
  }
  if((fout_cell_para=fopen(argv[optind+5], "w")) == NULL){
    printf("cannot open Output_file2 (parameters related to cell)\n");
  }
  if((fout_ll=fopen(argv[optind+6], "w")) == NULL){
    printf("cannot open Output_file3 (log-likelihood)\n");
  }

  int gene_num = atoi(argv[optind+7]);
  int cell_num = atoi(argv[optind+8]);

  Continuous_OU_process OU(gene_num, cell_num, K, max_ite1, max_ite2, alpha_min, alpha_max, t_min, t_max, sigma_squared_min, thresh);

  OU.Set_expression(fin_expression);
  if(OU.Set_initial_parameter(fin_init) == 1){
    return 1;
  }
  if(OU.Set_optimized_parameter(fin_opt_gene_para, fin_opt_cell_para) == 1){
    return 1;
  }

  OU.EM();

  OU.Print_cell_parameter(fout_cell_para);
  OU.Print_gene_parameter(fout_gene_para);
  OU.Print_ll(fout_ll);

  return 0;
}

int main_correlation(int argc, char* argv[]){
  int K = 1;

  int ch;
  extern char *optarg;
  extern int optind, opterr;
  while((ch = getopt(argc, argv, "k:")) != -1){
    switch(ch){
    case 'k':
      K = atoi(optarg);
      break;
    case ':':
      printf("invalid option\n");
      return 1;
    case '?':
      printf("invalid option\n");
      return 1;
    }
  }

  FILE *fin_expression, *fin_init, *fin_gene_para, *fin_cell_para, *fout_nexp, *fout_cor;
  if((fin_expression=fopen(argv[1], "r")) == NULL){
    printf("cannot open expression data\n");
    return 1;
  }
  if((fin_init=fopen(argv[2], "r")) == NULL){
    printf("cannot open initial distribution data\n");
    return 1;
  }
  if((fin_gene_para=fopen(argv[3], "r")) == NULL){
    printf("cannot open gene and lineage parameters data\n");
    return 1;
  }
  if((fin_cell_para=fopen(argv[4], "r")) == NULL){
    printf("cannot open cell parameters data\n");
    return 1;
  }
  if((fout_nexp=fopen(argv[5], "w")) == NULL){
    printf("cannot open Output_file1 (normalized expression)\n");
    return 1;
  }
  if((fout_cor=fopen(argv[6], "w")) == NULL){
    printf("cannot open Output_file2 (correlation matrix)\n");
    return 1;
  }

  int gene_num = atoi(argv[7]);
  int cell_num = atoi(argv[8]);

  Continuous_OU_process OU(gene_num, cell_num, K);

  OU.Set_expression(fin_expression);
  if(OU.Set_initial_parameter(fin_init) == 1){
    return 1;
  }
  if(OU.Set_optimized_parameter(fin_gene_para, fin_cell_para) == 1){
    return 1;
  }

  //calc normalized expression
  OU.Print_correlation(fout_nexp, fout_cor);

  return 0;
}

int main_sp(int argc, char* argv[]){
  FILE *fin_expression, *fin_initial, *fout_pseudo_time, *fout_mst, *fout_pca;
  if((fin_expression=fopen(argv[1], "r")) == NULL){
    printf("cannot open expression data\n");
    return 1;
  }
  if((fin_initial=fopen(argv[2], "r")) == NULL){
    printf("cannot open initial distribution data\n");
    return 1;
  }
  if((fout_pseudo_time=fopen(argv[3], "w")) == NULL){
    printf("cannot open Output_file1 (pseudo-time)\n");
    return 1;
  }
  if((fout_pca=fopen(argv[4], "w")) == NULL){
    printf("cannot open Output_file2 (PCA)\n");
    return 1;
  }

  int gene_num = atoi(argv[5]);
  int cell_num = atoi(argv[6]);
  int dim = atoi(argv[7]);

  Pseudo_Time PT(gene_num, cell_num, dim);
  PT.Set_expression(fin_expression);
  PT.Set_initial_parameter(fin_initial);

  //Prim
  PT.Prim(fout_pseudo_time, fout_pca);

  return 0;
}

#include <Rcpp.h>
// [[Rcpp::export]]
int main_wrap(Rcpp::CharacterVector args, bool verbose, int method) {
  if (!verbose) {
    fclose(stdout);
    fclose(stderr);
  }

  int argc = args.length();

  char **argv = (char**) malloc(argc * sizeof(char*));
  for (int i = 0; i < argc; i++) {
    std::string buf = Rcpp::as<std::string>(args(i));
    char *cstr = new char[buf.length() + 1];
    strcpy(cstr, buf.c_str());
    argv[i] = cstr;
  }

  if (method == 1) {
    return main_scoup(argc, argv);
  } else if (method == 2) {
    return main_scoup_resume(argc, argv);
  } else if (method == 3) {
    return main_correlation(argc, argv);
  } else if (method == 4) {
    return main_sp(argc, argv);
  } else {
    return 1;
  }
}
