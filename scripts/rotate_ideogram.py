import argparse
from PyPDF2 import PdfReader
from PyPDF2 import PdfFileWriter

def read_and_write_image(image, out):
    reader = PdfReader(image)
    page = reader.getPage(0)
    page.rotateClockwise(90)
    pdf_writer = PdfFileWriter()
    pdf_writer.addPage(page)
    pdf_out = open(out, 'wb')
    pdf_writer.write(pdf_out)
    pdf_out.close()   

def rotate_image(image1, image2, image3, image4, out1, out2, out3, out4):
    """
    """
    for image,out in zip([image1, image2, image3, image4], [out1, out2, out3, out4]):
        read_and_write_image(image, out)

def main():
    parser = argparse.ArgumentParser(
        """Function rotates pdf karyoplot images"""
    )
    parser.add_argument('--image1', type=str, help="relative path to VK ISCN page 1 of karyoplot")
    parser.add_argument('--image2', type=str, help="relative path to VK ISCN page 2 of karyoplot")
    parser.add_argument('--image3', type=str, help="relative path to VK ISCN page 3 of karyoplot")
    parser.add_argument('--image4', type=str, help="relative path to VK ISCN page 4 of karyoplot")
    parser.add_argument('--out1', type=str, help="relative path to rotated image")
    parser.add_argument('--out2', type=str, help="relative path to rotated image")
    parser.add_argument('--out3', type=str, help="relative path to rotated image")
    parser.add_argument('--out4', type=str, help="relative path to rotated image")

    args = parser.parse_args()
    print(args)
    image1 = args.image1
    image2 = args.image2
    image3 = args.image3
    image4 = args.image4
    out1 = args.out1
    out2 = args.out2
    out3 = args.out3
    out4 = args.out4
    rotate_image(image1=image1, image2=image2, image3=image3, image4=image4, out1=out1, out2=out2, out3=out3, out4=out4)

if __name__ == "__main__":
    main()