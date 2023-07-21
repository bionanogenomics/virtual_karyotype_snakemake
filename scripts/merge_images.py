import argparse
from PyPDF2 import PdfFileReader, PdfFileWriter
from PyPDF2 import PageObject


def merge_image(image1, image2, image3, image4, out1):
    """
    """
    page_1 = PdfFileReader(open(image1, 'rb')).getPage(0)
    page_2 = PdfFileReader(open(image2, 'rb')).getPage(0)
    page_3 = PdfFileReader(open(image3, 'rb')).getPage(0)
    page_4 = PdfFileReader(open(image4, 'rb')).getPage(0)

    translated_page = PageObject.create_blank_page(None, page_1.mediabox.width * 4, page_1.mediabox.height)
    translated_page.mergePage(page_1)
    translated_page.mergeScaledTranslatedPage(page_2, 1, page_2.mediabox.width, 0)
    translated_page.mergeScaledTranslatedPage(page_3, 1, page_3.mediabox.width * 2, 0)
    translated_page.mergeScaledTranslatedPage(page_4, 1, page_4.mediabox.width * 3, 0)
    translated_page.rotateClockwise(90)

    writer = PdfFileWriter()
    writer.addPage(translated_page)

    with open(out1, 'wb') as f:
        writer.write(f)

def main():
    parser = argparse.ArgumentParser(
        """Function rotates pdf karyoplot images"""
    )
    parser.add_argument('--image1', type=str, help="relative path to VK ISCN page 1 of karyoplot")
    parser.add_argument('--image2', type=str, help="relative path to VK ISCN page 2 of karyoplot")
    parser.add_argument('--image3', type=str, help="relative path to VK ISCN page 3 of karyoplot")
    parser.add_argument('--image4', type=str, help="relative path to VK ISCN page 4 of karyoplot")
    parser.add_argument('--out1', type=str, help="relative path to rotated image")


    args = parser.parse_args()
    print(args)
    image1 = args.image1
    image2 = args.image2
    image3 = args.image3
    image4 = args.image4
    out1 = args.out1

    merge_image(image1=image1, image2=image2, image3=image3, image4=image4, out1=out1)

if __name__ == "__main__":
    main()