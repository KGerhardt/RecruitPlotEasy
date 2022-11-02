#from recruit_plot_easy_2_database import run_build
#from recruit_plot_easy_2_plot import run_plot

import sys

def main():
	if len(sys.argv) < 2:
		sys.exit("Please choose a module: either 'build' or 'plot'")
	else:
		action = sys.argv[1]
		if action not in ["build", "plot"]:
			print("I could not recognize your module:", action+".", "Please try again.")
			sys.exit()
		
		if action == "build":
			from recruit_plot_easy_2_database import run_build
			run_build()
			
		if action == "plot":
			from recruit_plot_easy_2_plot import run_plot
			run_plot()
			
			
if __name__ == "__main__":
	main()