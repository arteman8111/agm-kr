export class Renderer {
    constructor(appSelector) {
        this.app = document.querySelector(appSelector);
    }

    createTitle(titleText) {
        return `<h2 style="margin: 25px 0 0 0; font-size: 24px; padding:10px 5px;">${titleText}</h2>`;
    }

    createTable(data, labels, startIndex) {
        const table = document.createElement('table');
        table.className = "data";

        const header = document.createElement('tr');
        header.innerHTML = `<td class="data__element"></td>`;
        
        if (startIndex !== 2) {
            header.innerHTML += data[0].map((_, i) => `<td class="data__element">${startIndex + i}</td>`).join('');
        }
        
        table.appendChild(header);

        data.forEach((row, i) => {
            const tr = document.createElement('tr');
            tr.innerHTML = `<td class="data__element">${labels[i]}</td>` +
                           row.map(element => `<td class="data__element">${math.format(element, { notation: 'exponential', precision: 5 })}</td>`).join('');
            table.appendChild(tr);
        });

        return table;
    }

    render(data, labels, startIndex, titleText = 'Заголовок') {
        const box = document.createElement('div');
        box.innerHTML = this.createTitle(titleText);
        box.appendChild(this.createTable(data, labels, startIndex));
        this.app.appendChild(box);
    }
}

export class InfoDisplay {
    constructor(appSelector) {
        this.app = document.querySelector(appSelector);
    }

    createInfoItem(term, value) {
        return `
            <div class="id__box">
                <dt class="id__item">${term}:</dt>
                <dd class="id__item">${value}</dd>
            </div>
        `;
    }

    render(group, N, h, M_inf, alfa, T_st, Re_k, c, b) {
        const infoHtml = `
            <section class="id">
                <h2 class="id__title">Мои ИД:</h2>
                <dl class="id__list">
                    ${this.createInfoItem('Группа', `СМ3-7${group}`)}
                    ${this.createInfoItem('Вариант N', N)}
                    ${this.createInfoItem('h, [м]', h)}
                    ${this.createInfoItem('M, [-]', M_inf)}
                    ${this.createInfoItem('α, [град]', math.round(alfa * 180 / 3.14))}
                    ${this.createInfoItem('T, [K]', T_st)}
                    ${this.createInfoItem('Reкр, [-]', Re_k)}
                    ${this.createInfoItem('с, [м]', c)}
                    ${this.createInfoItem('b, [м]', b)}
                </dl>
            </section>
        `;
        this.app.insertAdjacentHTML('afterbegin', infoHtml);
    }
}