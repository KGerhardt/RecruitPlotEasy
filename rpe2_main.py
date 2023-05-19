import sys

def main():
	if len(sys.argv) < 2:
		sys.exit("Please choose a module: either 'build' or 'plot' or 'describe'")
	else:
		action = sys.argv[1]
		if action not in ["build", "plot", "describe"]:
			print("I could not recognize your module:", action+".", "Please try again.")
			sys.exit()
		
		if action == "build":
			try:
				from .recruit_plot_easy_2_database import run_build
			except:
				from recruit_plot_easy_2_database import run_build
			run_build()
			
		if action == "plot":
			try:
				from .recruit_plot_easy_2_plot import run_plot
			except:
				from recruit_plot_easy_2_plot import run_plot
			run_plot()
			
			
		if action == "describe":
			try:
				from .recruit_plot_easy_2_describe import run_descriptor
			except:
				from recruit_plot_easy_2_describe import run_descriptor
			
			run_descriptor()
			
if __name__ == "__main__":
	main()