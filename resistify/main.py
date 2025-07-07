from resistify.cli import parse_args
from resistify._loguru import logger
from resistify.__version__ import __version__


def main():
    args = parse_args()

    if args.debug:
        logger.update_level("DEBUG")

    logger.info(f"Welcome to Resistify version {__version__}!")

    # Using lazy imports as it's pretty slow otherwise
    match args.command:
        case "nlr":
            from resistify.nlr import nlr

            nlr(args)
        case "prr":
            from resistify.prr import prr

            prr(args)
        case "download_models":
            from resistify.download_models import download_models

            download_models(args)
        case "draw":
            from resistify.draw import draw

            draw(args)

    logger.info("Thank you for using Resistify!")
    logger.info("If you used Resistify in your research, please cite the following:")
    logger.info(" - Resistify: https://doi.org/10.1177/11779322241308944")
    logger.info(" - NLRexpress: https://doi.org/10.3389/fpls.2022.975888")
    if args.command == "prr":
        logger.info(" - TMbed: https://doi.org/10.1186/s12859-022-04873-x")
    elif getattr(args, "coconat", False):
        logger.info(" - CoCoNat: https://doi.org/10.1093/bioinformatics/btad495")


if __name__ == "__main__":
    main()
